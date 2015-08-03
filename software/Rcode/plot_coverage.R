## ----- plot_coverage
## <<plot_coverage.R>>

## Plot the coverage for each sample and a given gene
## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

## path_to_display_file should be the path to a text file with five columns: 
## the first columns give the genes that are to be plotted, the second gives 
## the features that should be colored in green for each gene in gtf track 1, 
## the third gives the features that should be colored in red in gtf track 1, 
## and the fourth and fifth give the same for gtf track 2. The text file should 
## not have a header line. Multiple features to color should be separated by
## a comma.
print(path_to_tophat_output)
print(path_to_gtf_file)
print(path_to_second_gtf_file)
print(path_to_display_file)
print(highlight_type)
print(conditions)
print(pattern)
print(output_filename)
print(plot_height)
print(bounding_box)
print(show_name)
print(ext_before)
print(ext_after)

## Read the file defining which genes to display and which feature to highlight
displays <- read.delim(path_to_display_file, header = FALSE, 
                       colClasses = "character", as.is = TRUE)
genes_to_display <- displays[, 1]

highlight_green_list_gtf1 <- as.list(displays[, 2])
names(highlight_green_list_gtf1) <- genes_to_display
highlight_green_list_gtf1 <- lapply(highlight_green_list_gtf1, function(w) {
  paste0(gsub(" ", "", strsplit(w, ",")[[1]]), "_1")
})
highlight_red_list_gtf1 <- as.list(displays[, 3])
names(highlight_red_list_gtf1) <- genes_to_display
highlight_red_list_gtf1 <- lapply(highlight_red_list_gtf1, function(w) {
  paste0(gsub(" ", "", strsplit(w, ",")[[1]]), "_1")
})

highlight_green_list_gtf2 <- as.list(displays[, 4])
names(highlight_green_list_gtf2) <- genes_to_display
highlight_green_list_gtf2 <- lapply(highlight_green_list_gtf2, function(w) {
  paste0(gsub(" ", "", strsplit(w, ",")[[1]]), "_2")
})
highlight_red_list_gtf2 <- as.list(displays[, 5])
names(highlight_red_list_gtf2) <- genes_to_display
highlight_red_list_gtf2 <- lapply(highlight_red_list_gtf2, function(w) {
  paste0(gsub(" ", "", strsplit(w, ",")[[1]]), "_2")
})

## Help function to list subdirectories of a given directory
## From http://stackoverflow.com/questions/4749783/how-to-obtain-a-list-of-directories-within-a-directory-like-list-files-but-i
list.dirs <- function(path = ".", pattern = NULL, all.dirs = FALSE,
                      full.names = FALSE, ignore.case = FALSE) {
  all <- list.files(path, pattern, all.dirs,
                    full.names = TRUE, recursive = FALSE, ignore.case)
  dirs <- all[file.info(all)$isdir]
  if(isTRUE(full.names))
    return(dirs)
  else
    return(basename(dirs))
}

## List subdirectories in tophat directory
dir_names <- list.dirs(path_to_tophat_output, pattern = pattern, full.names = FALSE)
sample_names <- gsub("sample", "", dir_names)
bw_files <- paste0(path_to_tophat_output, "/", dir_names, "/accepted_hits.bw")

keep <- which(file.exists(bw_files))
dir_names <- dir_names[keep]
sample_names <- sample_names[keep]
bw_files <- bw_files[keep]

print(bw_files)
print(sample_names)

library(Gviz)
options(ucscChromosomeNames = FALSE)

library(rtracklayer)
genemodels <- import(path_to_gtf_file)
if (length(grep("\\.gtf$", path_to_gtf_file)) > 0) {
  idx <- match(c("transcript_id", "gene_id", "exon_id"), colnames(mcols(genemodels)))
  colnames(mcols(genemodels))[idx] <- c("transcript", "gene", "exon")
  mcols(genemodels)$symbol <- mcols(genemodels)$transcript
  genemodels <- subset(genemodels, type == "exon")
} else if (length(grep("\\.gff$", path_to_gtf_file)) > 0) {
  genemodels <- subset(genemodels, type == "exonic_part")
  mcols(genemodels)$exon <- sapply(as.character(mcols(genemodels)$group), 
                                   function(s) {
                                     m <- strsplit(s, ";")[[1]]
                                     m <- m[grep("exonic_part_number", m)]
                                     gsub("\"", "", 
                                          gsub("^\"", "\"E", 
                                               gsub(" exonic_part_number ", "", m)))
                                   })
  mcols(genemodels)$transcript <- sapply(as.character(mcols(genemodels)$group), 
                                         function(s) {
                                           m <- strsplit(s, ";")[[1]]
                                           m <- m[grep("gene_id", m)]
                                           gsub("\"", "", 
                                                gsub(" gene_id ", "", m))
                                         })
  mcols(genemodels)$gene <- mcols(genemodels)$transcript
  mcols(genemodels)$symbol <- mcols(genemodels)$transcript
} else {
  stop("path_to_gtf_file must lead to either a gtf or gff file")
}

if (!is.null(path_to_second_gtf_file)) {
  genemodels2 <- import(path_to_second_gtf_file)
  if (length(grep("\\.gtf$", path_to_second_gtf_file)) > 0) {
    idx <- match(c("transcript_id", "gene_id", "exon_id"), colnames(mcols(genemodels2)))
    colnames(mcols(genemodels2))[idx] <- c("transcript", "gene", "exon")
    mcols(genemodels2)$symbol <- mcols(genemodels2)$transcript
    genemodels2 <- subset(genemodels2, type == "exon")
  } else if (length(grep("\\.gff$", path_to_second_gtf_file)) > 0) {
    genemodels2 <- subset(genemodels2, type == "exonic_part")
    mcols(genemodels2)$exon <- sapply(as.character(mcols(genemodels2)$group), 
                                      function(s) {
                                        m <- strsplit(s, ";")[[1]]
                                        m <- m[grep("exonic_part_number", m)]
                                        gsub("\"", "", 
                                             gsub("^\"", "\"E", 
                                                  gsub(" exonic_part_number ", "", m)))
                                      })
    mcols(genemodels2)$transcript <- sapply(as.character(mcols(genemodels2)$group), 
                                            function(s) {
                                              m <- strsplit(s, ";")[[1]]
                                              m <- m[grep("gene_id", m)]
                                              gsub("\"", "", 
                                                   gsub(" gene_id ", "", m))
                                            })
    mcols(genemodels2)$gene <- mcols(genemodels2)$transcript
    mcols(genemodels2)$symbol <- mcols(genemodels2)$transcript
  } else {
    stop("path_to_second_gtf_file must lead to either a gtf or gff file")
  }
}

pdf(output_filename, height = plot_height, width = 7)
for (gene_to_display in genes_to_display) {
  if (isTRUE(bounding_box)) col <- "black"
  else col <- NULL
  
  ## Prepare first gene model track
  gm <- subset(genemodels, gene == gene_to_display)
  show_chr <- unique(seqnames(gm))[1]
  min_coord <- min(start(gm)) - ext_before*(max(end(gm)) - min(start(gm)))
  max_coord <- max(end(gm)) + ext_after*(max(end(gm)) - min(start(gm)))
  grtr <- GeneRegionTrack(gm, showId = show_name, name = gene_to_display,
                          col = col, fill = "orange")
  group(grtr) <- mcols(gm)$transcript
  
  if (highlight_type == "transcript") {
    feature(grtr) <- paste0(mcols(gm)$transcript, "_1")
  } else if (highlight_type == "exon") {
    feature(grtr) <- paste0(mcols(gm)$exon, "_1")
  }
  
  ## If a second gtf file is provided, prepare second gene model track
  if (!is.null(path_to_second_gtf_file)) {
    gm2 <- subset(genemodels2, gene == gene_to_display)
    grtr2 <- GeneRegionTrack(gm2, showId = show_name, name = gene_to_display,
                             col = col, fill = "orange")
    group(grtr2) <- mcols(gm2)$transcript
    if (highlight_type == "transcript") {
      feature(grtr2) <- paste0(mcols(gm2)$transcript, "_2")
    } else if (highlight_type == "exon") {
      feature(grtr2) <- paste0(mcols(gm2)$exon, "_2")
    }
  }
  
  gtr <- GenomeAxisTrack()
  
  twocols <- c(rgb(182, 109, 255, maxColorValue = 255), 
               rgb(219, 209, 0, maxColorValue = 255))
  
  multiTracks <- lapply(1:length(bw_files), function(i) {
    assign(paste0("dtr", i), 
           DataTrack(range = bw_files[i],
                     type = "histogram",
                     name = sample_names[i],
                     chromosome = unique(seqnames(gm)),
                     fill = twocols[(as.numeric(as.factor(conditions)))[i]],
                     col = twocols[(as.numeric(as.factor(conditions)))[i]],
                     col.histogram = twocols[(as.numeric(as.factor(conditions)))[i]],
                     fill.histogram = twocols[(as.numeric(as.factor(conditions)))[i]]))
  })
  
  ## Build up the coloring
  myoptions <- ""
  if (!is.null(highlight_red_list_gtf1[[gene_to_display]])) {
    for (hr in highlight_red_list_gtf1[[gene_to_display]]) {
      myoptions <- paste0(myoptions, "'", hr, "' = 'red', ")
    }
  }
  if (!is.null(highlight_green_list_gtf1[[gene_to_display]])) {
    for (hg in highlight_green_list_gtf1[[gene_to_display]]) {
      myoptions <- paste0(myoptions, "'", hg, "' = 'green', ")
    }
  }
  
  if (!is.null(highlight_red_list_gtf2[[gene_to_display]])) {
    for (hr in highlight_red_list_gtf2[[gene_to_display]]) {
      myoptions <- paste0(myoptions, "'", hr, "' = 'red', ")
    }
  }
  if (!is.null(highlight_green_list_gtf2[[gene_to_display]])) {
    for (hg in highlight_green_list_gtf2[[gene_to_display]]) {
      myoptions <- paste0(myoptions, "'", hg, "' = 'green', ")
    }
  }
  myoptions <- gsub(", $", "", myoptions)
  #myoptions <- paste(highlight_red, "= 'red',", highlight_green, "= 'green'")
  
  tracks <- c(multiTracks, gtr, grtr)
  if (!is.null(path_to_second_gtf_file)) {
    tracks <- c(tracks, grtr2)
  }
  
  eval(parse(text = paste("plotTracks(tracks, 
                                   main = gene_to_display, 
                                   from = min_coord, to = max_coord, 
                                   min.width = 0, min.distance = 0, collapse = FALSE,", 
                          myoptions, ")")))
  
#   eval(parse(text = paste("plotTracks(tracks, 
#                                    col = NULL, main = gene_to_display, 
#                                    from = min_coord, to = max_coord, 
#                                    min.width = 0, min.distance = 0, collapse = FALSE,", 
#                           myoptions, ")")))
  
  # plotTracks(c(multiTracks, gtr, grtr), 
  #            col = NULL, from = min_coord, to = max_coord, 
  #            min.width = 0, min.distance = 0, collapse = FALSE)
}
dev.off()




