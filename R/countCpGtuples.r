# Peter Hickey
# 17/04/2012
# Count the number CpG-k-tuples in a genome within a maximum distance, where the maximum distance corresponds to the maximum read (fragment) length in single-end (paired-end) sequencing

library(rtracklayer)
x <- import.bed("hg19_CpGs.bed.gz") # Found in home area on unix boxes
library(GenomicRanges)
xx <- as(Class="GRanges", object = x)
rm(x)
gc()

# Returns a a count of the number of CpG-k-tuples within maxlength distance, e.g. countCpGTuples(2, xx, 100) is the number of CpG-pairs within 100bp distance
# k is the size of the tuple to be considered
# x is the list of CpGs in the genome as a GRanges object
# maxlength is the maximum separation of CpGs to be considered
countCpGTuples <- function(k, x, maxlength){
	d <- diff(start(x), lag = (k-1))
	n <- sum(d > 0 & d <= maxlength)
	return(n)
}

library(plyr)
require(doMC)
registerDoMC(cores = 15)
# Maximum separation of CpGs; defined by the read- or fragment-length, as appropriate.
maxlength <- 500
m <- maxlength/2
l <- 2:m
n <- aaply(.data = l, .fun = countCpGTuples, .margin = 1, x = xx, maxlength = maxlength, .parallel = TRUE)

plot(n)
plot(n/length(xx) * 100, xlab = "x", ylab = "%", main = paste("% of CpGs with x CpGs within ", maxlength, "bp", sep = ""), type = "h") 
plot(n/length(xx) * 100, xlab = "x", ylab = "%", main = paste("% of CpGs with x CpGs within ", maxlength, "bp", sep = ""), xlim = c(1, 50), type = "h") # The % of CpGs with x CpGs within maxlength (bp)

# Pretty plots with ggplot2
library(ggplot2)