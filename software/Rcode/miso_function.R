## ----- miso_function
## <<miso_function.R>>

## Summarise and merge the MISO output for all samples
miso <- function(miso_output_dir, counts_output_dir) {
  ## Sort genes with respect to gene ID
  allgenes <- dir(miso_output_dir, "miso$", recursive = TRUE, full.names = TRUE)
  gn <- gsub(".miso", "", basename(allgenes))
  
  samplename <- sapply(allgenes, function(u) {
    a <- strsplit(u, "/")[[1]]
    a[grep("sample", a)]
  })
  genename <- sapply(allgenes, function(u) {
    a <- strsplit(u, "/")[[1]]
    a[grep("\\.miso$", a)]
  })
  names(allgenes) <- paste0(samplename, ":", genename)
  
  ## For each .miso file, create a named vector with the
  ## event (names) and the count (numeric)
  ## First check that the isoforms are in the same order for all samples
  A <- lapply(allgenes, function(u) {
    read.table(u, comment.char = "", nrow = 1, sep = "\t", 
               stringsAsFactors = FALSE)$V1
  })
  As <- split(A, gn)
  stopifnot(all(sapply(As, function(w) length(unique(w))) == 1))

  g <- lapply(allgenes, function(u) {
    rt <- read.table(u, comment.char = "", nrow = 1, sep = "\t", 
                     stringsAsFactors = FALSE)$V8
    gsub("counts=\\(", "", rt)
  })
  h <- lapply(g, function(u) strsplit(u, ",\\(")[[1]])
  j <- lapply(h, function(u) {
    z <- unlist(strsplit(u, "\\):"))
    n <- length(z)
    x <- as.numeric(z[seq(2, n, by = 2)])
    names(x) <- z[seq(1, n, by = 2)]
    x
  })
  
  ## Break into a list for each gene
  js <- split(j, gn)
  jsl <- sapply(js, length)
  
  js2 <- lapply(js, function(w) {
    if (length(w) != max(jsl)) {
      pres <- sapply(strsplit(names(w), ":"), .subset, 1)
      miss <- setdiff(paste0("sample", 1:max(jsl)), pres)
      gen <- strsplit(names(w)[1], ":")[[1]][2]
      for (i in 1:length(miss)) {
        w[[paste0(miss[i], ":", gen)]] <- 0
      }
    }
    w
  })
  
  ## Function to go from list to table
  getTableFromList <- function(u) {
    cn <- sapply(names(u), function(v) strsplit(v,":")[[1]][1], USE.NAMES = FALSE)
    rn <- unique(unlist(lapply(u, names)))
    x <- matrix(0, nrow = length(rn), ncol = length(cn), dimnames = list(rn, cn))
    for(i in 1:length(u))
      x[names(u[[i]]), cn[i]] <- u[[i]]
    x
  }
  
  ks <- lapply(js2, getTableFromList)
  ks <- lapply(ks, function(w) w[, order(colnames(w)), drop = FALSE])
  
  n <- sapply(ks, nrow)
  rn <- rep(names(ks), n)
  
  ## Create one big table
  MISO_count <- do.call("rbind", ks)
  rownames(MISO_count) <- paste0(rn, ":", rownames(MISO_count))
  MISO_count <- data.frame(event = rownames(MISO_count), MISO_count)
  
  snames <- sapply(colnames(MISO_count), function(u) gsub("sample", "", u))
  for(i in 2:ncol(MISO_count)){
    final.counts <- data.frame(MISO_count[, 1], MISO_count[, i])  
    write.table(final.counts, paste0(counts_output_dir, "/miso", 
                                     snames[i], ".txt"),
                col.names = FALSE, row.names = FALSE, sep = "\t",
                quote = FALSE) 
  } 
}

## Summarise and merge the MISO output for all samples.
## Transcript level counts
miso_transcript <- function(miso_output_dir, counts_output_dir) {
  ## Sort genes with respect to gene ID
  allgenes <- dir(miso_output_dir, "miso$", recursive = TRUE, full.names = TRUE)
  gn <- gsub(".miso", "", basename(allgenes))
  
  samplename <- sapply(allgenes, function(u) {
    a <- strsplit(u, "/")[[1]]
    a[grep("sample", a)]
  })
  genename <- sapply(allgenes, function(u) {
    a <- strsplit(u, "/")[[1]]
    a[grep("\\.miso$", a)]
  })
  names(allgenes) <- paste0(samplename, ":", genename)
  
  ## For each .miso file, create a named vector with the
  ## event (names) and the count (numeric)
  ## First check that the isoforms are in the same order for all samples
  A <- lapply(allgenes, function(u) {
    read.table(u, comment.char = "", nrow = 1, sep = "\t", 
               stringsAsFactors = FALSE)$V1
  })
  As <- split(A, gn)
  stopifnot(all(sapply(As, function(w) length(unique(w))) == 1))
  
  g <- lapply(allgenes, function(u) {
    rt <- read.table(u, comment.char = "", nrow = 1, sep = "\t", 
                     stringsAsFactors = FALSE)$V9
    gsub("assigned_counts=", "", rt)
  })
  h <- lapply(g, function(u) strsplit(u, ",")[[1]])
  j <- lapply(h, function(u) {
    z <- unlist(strsplit(u, ":"))
    n <- length(z)
    x <- as.numeric(z[seq(2, n, by = 2)])
    names(x) <- z[seq(1, n, by = 2)]
    x
  })
  
  ## Break into a list for each gene
  js <- split(j, gn)
  jsl <- sapply(js, length)
  
  js2 <- lapply(js, function(w) {
    if (length(w) != max(jsl)) {
      pres <- sapply(strsplit(names(w), ":"), .subset, 1)
      miss <- setdiff(paste0("sample", 1:max(jsl)), pres)
      gen <- strsplit(names(w)[1], ":")[[1]][2]
      for (i in 1:length(miss)) {
        w[[paste0(miss[i], ":", gen)]] <- 0
      }
    }
    w
  })
  
  ## Function to go from list to table
  getTableFromList <- function(u) {
    cn <- sapply(names(u), function(v) strsplit(v,":")[[1]][1], USE.NAMES = FALSE)
    rn <- unique(unlist(lapply(u, names)))
    x <- matrix(0, nrow = length(rn), ncol = length(cn), dimnames = list(rn, cn))
    for(i in 1:length(u))
      x[names(u[[i]]), cn[i]] <- u[[i]]
    x
  }
  
  ks <- lapply(js2, getTableFromList)
  ks <- lapply(ks, function(w) w[, order(colnames(w)), drop = FALSE])
  
  n <- sapply(ks, nrow)
  rn <- rep(names(ks), n)
  
  ## Create one big table
  MISO_count <- do.call("rbind", ks)
  rownames(MISO_count) <- paste0(rn, ":", rownames(MISO_count))
  MISO_count <- data.frame(event = rownames(MISO_count), MISO_count)
  
  snames <- sapply(colnames(MISO_count), function(u) gsub("sample", "", u))
  for(i in 2:ncol(MISO_count)){
    final.counts <- data.frame(MISO_count[, 1], MISO_count[, i])  
    write.table(final.counts, paste0(counts_output_dir, "/miso", 
                                     snames[i], ".txt"),
                col.names = FALSE, row.names = FALSE, sep = "\t",
                quote = FALSE) 
  } 
}
