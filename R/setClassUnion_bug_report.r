# Reproducible example of problem with setClassUnion("vectorOrNULL", c("vector", "NULL")) (defined in WFComethylation_class.r)

## This version works as intended
library(GenomicRanges)

out <- list(chr = rep('chr21', 10), 1:10, start = 1:10, end = 2:11)
showMethods('Rle')
gr <- GRanges(seqnames = out[['chr']], ranges = IRanges(start = out[['start']], end = out[['end']]))
gr
showMethods('Rle')
sessionInfo()

## But this version does not
## Firstly, start a fresh R session
library(GenomicRanges)
setClassUnion("vectorOrNULL", c("vector", "NULL")) ## This line is the culprit 
out <- list(chr = rep('chr21', 10), 1:10, start = 1:10, end = 2:11)
showMethods('Rle')
gr <- GRanges(seqnames = out[['chr']], ranges = IRanges(start = out[['start']], end = out[['end']]))
showMethods('Rle')