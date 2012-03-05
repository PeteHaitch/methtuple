## A modification of tabulate_comethylation.R to use DataFrame's instead of data.frame's
# A script to tabulate the methylation state of reads overlapping two proximal CpGs
# Peter Hickey (hickey@wehi.edu.au)
# 05/12/2011

library(Rsamtools)
library(BSgenome)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)
library(multicore)

# Function to read and process the bismark BAM file
# Assume reads have been aligned to the hg19 reference genome
# Requires BSgenome and BSgenome.Hsapiens.UCSC.hg19 to have already been loaded
readBismarkBam <- function(x, reference = NA){
  bam.param <- ScanBamParam(what = c("rname", "strand", "pos", "qwidth"), tag = "XM")
  bam <- scanBam(file = x, index = x, param = bam.param)
  bam[[1]]$rname <- as(bam[[1]]$rname, "character")
  bam[[1]]$strand <- as(bam[[1]]$strand, "character")
  gr <- GRanges(seqnames=bam[[1]]$rname, strand = bam[[1]]$strand, ranges=IRanges(start=bam[[1]]$pos, width = bam[[1]]$qwidth), XM = bam[[1]]$tag$XM)
  # Add reference genome metadata
  seqlevels(gr) <- seqlevels(Hsapiens)
  seqlengths(gr) <- seqlengths(Hsapiens)
  gr
}

### Functions for processing the XM tag
## A function for the Watson strand reads, the Crick strand reads and a "master" function for all reads
## I use the forward/reverse Watson/Crick definitions from http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3055211/
## Functions take in a GRanges instance containing the reads (x) and a RangedData instance containing the pair of CpGs (cpg.pair)
## The cpg.pair object has start(cpg.pair) = C in CpG_1 and end(cpg.pair) = G in CpG_2

# Extract methylation call at CpGs in cpg.pair for reads (x) aligned to the Watson strand
# Will break if length(x) == 0
watsonMethylCall <- function(x, cpg.pair){
    y1 <- substr(elementMetadata(x)$XM, start(cpg.pair) - start(x) + 1, start(cpg.pair) - start(x) + 1)
    y2 <- substr(elementMetadata(x)$XM, end(cpg.pair) - start(x), end(cpg.pair) - start(x)) 
    return(DataFrame(first = y1, second = y2))
}

# Extract methylation call at CpGs in cpg.pair for reads (x) aligned to the Crick strand
# Will break if length(x) == 0
# The XM tag for a Crick-strand read uses "left-most" coordinates
crickMethylCall <- function(x, cpg.pair){
    y1 <- substr(elementMetadata(x)$XM, start(cpg.pair) - start(x) + 2, start(cpg.pair) - start(x) + 2)
    y2 <- substr(elementMetadata(x)$XM, end(cpg.pair) - start(x) + 1, end(cpg.pair) - start(x) + 1) 
    return(DataFrame(first = y1, second = y2))
}

# getMethylCall extracts the methylation call for reads overlapping a particular CpG pair
# cpg.pair.idx = An index relating the particular CpG pair being processed back to the cpg.pairs object
# overlaps.list = An index of reads that overlap the particular CpG pair being processed
# reads = The full set of reads for the sample being processed
# cpg.pairs = The full set of pairs of CpGs
# Assumes reads align to + or - strand, i.e. strand=="*" is not allowed
# Returns the results as a DataFrame
getMethylCall <- function(cpg.pair.idx, overlaps.list, reads, cpg.pairs){
  cpg.pair.idx <- as.numeric(cpg.pair.idx)
  reads <- reads[as.vector(unlist(overlaps.list[cpg.pair.idx]))]
  cpg.pair <- cpg.pairs[cpg.pair.idx, ]
  x <- split(reads, strand(reads), drop = TRUE)
  # Process reads aligned to Watson strand (if there are any such reads)
  if(length(x$'+') > 0){
    watson <- watsonMethylCall(x$'+', cpg.pair)
  } else{
    watson <- vector()
  }
  # Process reads aligned to Crick strand (if there are any such reads)
  if(length(x$'-') > 0){
    crick <- crickMethylCall(x$'-', cpg.pair)
  } else{
    crick <- vector()
  }
  # Combine the results and return as a DataFrame 
  # This is a kludge that stores the index of the cpg.pair in the row names of the output
   output <- rbind(watson, crick)
   rownames(output)[1] <- cpg.pair.idx
   return(output)
}

## Using findOverlaps to do the subsetting (only required to be run once per sample!)
# For each read, choose the first or "leftmost" CpG pair it overlaps
firstCpGPair <- function(reads, cpg.pairs){
  overlaps <- findOverlaps(query = reads, subject = cpg.pairs)
  idx <- !duplicated(as.matrix(overlaps)[,1])
  as.matrix(overlaps)[idx,]
}

# For each read choose a random CpG pair that it overlaps
randomCpGPair <- function(reads, cpg.pairs, verbose = FALSE){
  overlaps <- findOverlaps(query = reads, subject = cpg.pairs)
  # olaps.table counts the number of CpG pair each read overlaps (read=query)
  olaps.table <- table(queryHits(overlaps))
  # n.olpas is the number of CpG pair each read overlaps
  n.olaps <- as.vector(olaps.table)
  # x_i is randomly selected from 1, 2,..., n.olaps_i (i.e. each entry in olaps.table)
  # i = 1, 2, ..., n.unique.reads
  x <- apply(X = olaps.table, FUN = sample, MARGIN = 1, size = 1)
  # which.pair realtes x back to olaps.table (as a row-index)
  which.pair <- cumsum(n.olaps) - n.olaps + x
  
  ## Verbose output
  if(verbose){
    data.frame(query = as.numeric(names(olaps.table)), n.olaps = n.olaps,
             x = x, which.pair = which.pair, subject = as.matrix(overlaps)[which.pair,2])
  } else{
  ## Concise output
  as.matrix(overlaps)[which.pair,]
  }
}

# A wrapper function for selecting the CpG pairs for each read (either "leftmost" or random selection)
chooseCpGPair <- function(reads, cpg.pairs, type){
  if(type == "first"){
    firstCpGPair(reads, cpg.pairs)
  } else if(type == "random"){
    randomCpGPair(reads, cpg.pairs)
  } else{
    print("Error: type argument must be \'first\' or \'random\'")
  }
}

# Converts a "findOverlaps-style" matrix (e.g. output from randomCpGPair()) to a list
# Each entry of the list is a CpG pair containing a vector of the (filtered) reads overlapping that pair
# This is a list of lists, i.e. for each chromosome ("List") there is an entry ("list") for each CpG pair with > 0 overlapping reads 
# The "list" stores which (filtered) reads overlap that particular pair of CpGs
# List = chr1
#   list = pos_1
#   list = pos_n_chr1
# List = chr2
#   list = pos_1
# etc.
olaps2List <- function(x){
  # if length(x) == 0 then there are no CpG pairs on this chromosome/contig
  if(length(x) == 0){
    return(numeric(0))
  }
  # If length(x) == 2 then there is a single pair of CpGs on this chromosome/contig
  if(length(x) == 2){
    x <- matrix(x, ncol = 2)
  } else{
    if(!is.null(dim(x))){
      # Must be sorted by subject
      if(is.unsorted(x[,2])){
        x <- x[order(x[,2], decreasing = FALSE), ]
      } else{
      }
    }
  }
  ## The conversion process
  # Use rle() to record which reads overlap which pairs
  subject.rle <- rle(x[, 2])
  n <- length(subject.rle$lengths)
  reads.per.CpG.pair <- vector("list", n)
  row.idx <- 1
  subject.rle.lengths <- subject.rle$lengths
  for(i in seq_len(n)){
    reads.per.CpG.pair[[i]] <- x[seq.int(row.idx, row.idx + subject.rle.lengths[i] - 1, 1), 1]
    # Update row.idx to point to the next subject to be processed
    row.idx <- row.idx + subject.rle.lengths[i]
  }
  # Add names to the CpG pairs
  names(reads.per.CpG.pair) <- unique(x[,2])
  reads.per.CpG.pair
}

makeTable <- function(chr, reads, cpg.pairs, input.for.makeTable){ 
  chr.data <- input.for.makeTable[[chr]]
  cpg.pairs <- cpg.pairs[[chr]]
  mclapply(X = names(chr.data), FUN = getMethylCall, overlaps.list = chr.data, reads = reads, cpg.pairs = cpg.pairs)
}

# Computes the comethylation table for a single chromosome
# reads = output of readBismarkBam()
# cpgs = a bed file of CpGs in the reference genome imported with rtracklayer::import()
# pair.choice = For reads that overlap multiple CpG pairs select whether the "leftmost" pair or a "random" pair is chosen
coMethTable <- function(reads, cpgs, pair.choice = "random"){
  # Create pairs of CpGs
  # Stored as a list with each element of the list corresponding to a chromosome/contig
  params.cpg.pairs <- RDApplyParams(rangedData = bed, applyFun = makeCpGPairs)
  genome.cpg.pairs <- rdapply(x = params.cpg.pairs)
  
  # Find all overlaps between the reads and the CpG pairs, 
  # then select one CpG pair per read (the choice of pair is controlled by the pair.choice argument)
  # Results are stored in a list by chromosome/contig
  overlaps.list <- lapply(genome.cpg.pairs, chooseCpGPair, reads = reads, type = pair.choice)
  
  # Create input for the makeTable() function
  # This is a list of lists, i.e. for each chromosome ("List"") thereis an entry ("list") for each positions
  # the "list" stores which reads overlap that particular pair of CpGs
  # List = chr1
  #   list = pos_1
  #   list = pos_n_chr1
  # List = chr2
  #   list = pos_1
  # etc.
  input.for.makeTable <- lapply(overlaps.list, olaps2List)
  
  # chr.names stores the names of all the chromosomes to be processed
  chr.names <- names(input.for.makeTable)
  
  # Create the tables for every CpG pair on every chromosome
  output <- lapply(chr.names, makeTable, reads = reads, cpg.pairs = genome.cpg.pairs, input.for.makeTable = input.for.makeTable)
  names(output) <- names(cpg.pairs)
}

# makeCpGPairs() transforms a RangedData object of CpGs to a RangedData object of neighbouring CpG pairs
# A makeCpGPairs() object stores the C of CpG_1 in the start() and the G of CpG_2 in the end() of the range
makeCpGPairs <- function(cpgs, type = "disjoint"){
  # Create "disjoint" CpG pairs, e.g. (1, 2), (3, 4), ..., (n-3, n-2), (n-1, n) 
  disjoint <- function(x){
    n <- nrow(x)
    if(n %% 2 == 1){
      n <- n - 1
    }
    RangedData(ranges = IRanges(start = start(x)[seq.int(1, n, 2)], 
                              end = start(x)[seq.int(2, n, 2)] + 1), space = space(x)[seq.int(1, n, 2)])      
    }
  # Create "overlapping" CpG pairs, e.g. (1, 2), (2, 3), ..., (n-2, n-1), (n-1, n)
  overlapping <- function(x){
    n <- nrow(x)
    RangedData(ranges = IRanges(start = start(x)[1:(n-1)], 
                              end = start(x)[2:n] + 1), space = space(x)[1:(n-1)])
  }
  
  if(type == "disjoint"){
      params.cpg.pairs <- RDApplyParams(rangedData = cpgs, applyFun = disjoint)
      rdapply(x = params.cpg.pairs)
  } else if(type == "overlapping"){
      params.cpg.pairs <- RDApplyParams(rangedData = cpgs, applyFun = overlapping)
      rdapply(x = params.cpg.pairs)    
  } else{
    print("Error: \'type\'' must be one of \'disjoint\' or \'overlapping\'.")
  }
}

## Computes the comethylation tables for a single sample
## reads = output of readBismarkBam()
## pair.choice = For reads that overlap multiple CpG pairs select whether the "leftmost" pair or a "random" pair is chosen
makeCoMethTable <- function(reads, cpg.pairs, pair.choice = "random"){
  # Find all overlaps between the reads and the CpG pairs, then select one CpG pair per read (the choice of pair is controlled by the pair.choice argument)
  # The result is a list of RangesMatchingList instances, where each element of the list corresponds to a chromosome/contig
  # The "query" in the RangesMatchList instances correspond to reads and the "subject" to CpG pairs
  overlaps.list <- lapply(cpg.pairs, chooseCpGPair, reads = reads, type = pair.choice)
  
  # Create input for the makeTable() function
  input.for.makeTable <- lapply(overlaps.list, olaps2List)
  
  # chr.names stores the names of all the chromosomes to be processed
  chr.names <- names(input.for.makeTable)
  
  # Create the tables for every CpG pair on every chromosome
  tab <- lapply(chr.names, makeTable, reads = reads, cpg.pairs = cpg.pairs, input.for.makeTable = input.for.makeTable)
  #names(tab) <- names(cpg.pairs)
  names(tab) <- chr.names
  
  ## Now I want to add those comethylation results as metadata to the corresponding CpG pair
  # First convert cpg.pairs from a RangedList to a GRanges object to allow insertion of metadata
  cpg.pairs <- mclapply(X = cpg.pairs, FUN = as, Class = "GRanges")
  # getID() is a function to extract the index for the cpg.pair (the ID) from the comethylation table
  getID <- function(x){
    unlist(lapply(x, function(x){as.numeric(rownames(x)[1])})) ### THIS WILL PROBABLY BREAK
  }
  
  # addComethData() adds as metadata the corresponding comethylation table to each CpG pair
  # chr.name = the name of the chromosome to be processed
  # cpg.pairs = the pairs of CpGs stored as a list of GRanges object
  # cometh.tab = the comethylation tables for each chromosome stored as a list
  addComethData <- function(chr.name, cpg.pairs, cometh.tab){
   x <- cpg.pairs[[chr.name]]
    y <- cometh.tab[[chr.name]]
    IDs <- getID(y)
    cometh.data <- vector("list", length(x)) # Needs to be a list of data.frames
    cometh.data[IDs] <- y
    x@elementMetadata$cometh.data <- cometh.data
    return(x)
   }
  #cpg.pairs.with.cometh <- unlist(lapply(names(cpg.pairs), addComethData, cpg.pairs = cpg.pairs, cometh.tab = tab))
  cpg.pairs.with.cometh <- unlist(lapply(chr.names, addComethData, cpg.pairs = cpg.pairs, cometh.tab = tab))
  names(cpg.pairs.with.cometh) <- chr.names
  return(cpg.pairs.with.cometh)
}

### QUESTIONS
# Should I only be using each read once? YES
# Should each CpG only contribute to one CpG pair? YES (for now)
# How does these two restrictions affect one another
# Can performance be improved by storing the reads and the CpG Pairs as IntervalTree objects?