# Peter Hickey
# 05/07/2012
# Aggregated-by-distance within-fragment comethylation plots, stratified by CGI-status and average methylation

#### TODOs ####
# Save plots as pdfs
# Save plots at "Powerpoint resolution"
# Add labels, titles, etc. to plots
# Decide whether it is better to stratify by beta or by gamma

#### Read in command line arguments ####
args <- commandArgs(TRUE)
# The name of the sample of interest
sample.name <- args[1]

#### Load libraries ####
library(stringr)
library(GenomicRanges)
library(plyr)
library(ggplot2)
library(parallel)
library(reshape2)
library(BSgenome)
library('BSgenome.Hsapiens.UCSC.hg18')

#### Load WF and AM files ####
load(paste(sample.name, '_comethylation.RData', sep = ''))
# Remove unrequired objexts
rm(args, CGI.shore, countInterveningCpGs, knownGene, knownGene.exons, lor, min.cov, n.cores, nonexonic.knownGene, outside.CGI, outside.knownGene, outside.knownGene.exons, path.to.gf, plotAggregatedLorByLag, plotLorBoxplotByLag, readWF, TES_2000, TES_500, tmp.gr, tmp.lor, TSS_2000, TSS_500, WF.gr)
gc()
AM <- read.table(str_c('../AM/', sample.name, '.am'), header = TRUE, stringsAsFactors = FALSE, )
AM.gr <- GRanges(seqnames = AM$chr, ranges = IRanges(start = AM$pos, end = AM$pos + 1), beta = AM$beta, gamma = AM$gamma)
rm(AM)
gc()

#### Function definitions ####
# Compute the log-odds ratio (lor) for a single CpG-pair or aggregate across CpG-pairs and then compute lor
# x is the output of readWF(), strand is which strand is analysed (combined, OT or OB), aggregate is whether to aggregate the counts across multiple CpG-pairs
lor <- function(x, strand = 'combined', aggregateByLag = FALSE, aggregateByNIC = FALSE){
  if(aggregateByLag & aggregateByNIC){
    stop('Only ONE of aggregateByLag or aggregateByNIC can be TRUE')
  }
  x.em <- elementMetadata(x)
  if(aggregateByLag == TRUE){
    agg.df <- data.frame(lag = width(x) - 1,
                         MM = x.em$MM,
                         MU = x.em$MU,
                         UM = x.em$UM,
                         UU = x.em$UU
    )
    agg.counts <- ddply(.data = agg.df, .variables = 'lag', .fun = function(x){colSums(x[, c('MM', 'MU', 'UM', 'UU')])}, .parallel = TRUE)                     
    MM <- agg.counts$MM
    MU <- agg.counts$MU
    UM <- agg.counts$UM
    UU <- agg.counts$UU
  } else{
    if(strand == 'combined'){
      MM <- x.em$MM
      MU <- x.em$MU
      UM <- x.em$UM
      UU <- x.em$UU
    }else if(strand == 'OT'){
      MM <- x.em$MM_OT
      MU <- x.em$MU_OT
      UM <- x.em$UM_OT
      UU <- x.em$UU_OT
    }else if(strand == 'OB'){
      MM <- x.em$MM_OB
      MU <- x.em$MU_OB
      UM <- x.em$UM_OB
      UU <- x.em$UU_OB
    }
  }
  # Continuity correct the counts
  MM <- MM + 0.5
  MU <- MU + 0.5
  UM <- UM + 0.5
  UU <- UU + 0.5
  # Compute lor and its asymptotic standard error (ASE)
  lor <- ifelse(MM == 0.5 & MU == 0.5 & UM == 0.5 & UU == 0.5, NA, log(MM * UU / (MU * UM)))
  ASE <- ifelse(MM == 0.5 & MU == 0.5 & UM == 0.5 & UU == 0.5, NA, sqrt(1/MM + 1/MU + 1/UM + 1/MU))
  n <- ifelse(MM == 0.5 & MU == 0.5 & UM == 0.5 & UU == 0.5, NA, MM + MU + UM + UU - 4 * 0.5) # Remove continuity correction
  if(aggregateByLag == TRUE){
    return.df <- data.frame(lag = agg.counts$lag, lor = lor, ASE = ASE, n = n)
  } else{
    return.df <- data.frame(lor = lor, ASE = ASE)
  }
  return(return.df)
}

# Count the number of CpGs in a 'w' sized window around a CpG-pair
# 'genome' is a BSgenome object, 'x' is a GRanges object, e.g. 'merged', 'w' is the window size
# 'x' must have well-defined 'seqlengths' attribute. Otherwise 'resize' function cannot correctly trim regions with negative start values or end values > length of the chromosome
countCpGsInWindow <- function(x, genome, w){
  regions <- resize(x, fix = 'center', width = w)
  CpG.seq <- getSeq(genome, names = regions, as.character = FALSE)
  n.CpGs <- vcountPattern(pattern = "CG", subject = CpG.seq) 
}

# Compute the average methylation in a 'w' sized window around a CpG-pair
# 'WF.gr' is the WF.gr object that contains all CpG-pairs in the sample, 'w' is the window size, 'AM.gr' is the AM.gr object that contains the beta-values for all CpGs in the sample
# 'WF.gr' must have well-defined 'seqlengths' attribute. Otherwise 'resize' function cannot correctly trim regions with negative start values or end values > length of the chromosome
# Code based on https://stat.ethz.ch/pipermail/bioc-sig-sequencing/2010-June/001325.html
averageMethylationInWindow <- function(x, AM.gr, w){
  regions <- resize(x, fix = 'center', width = w)
  means <- data.frame(beta_w = rep(NA, length(regions)), gamma_w = rep(NA, length(regions)))
  row.index <- 1
  for(chrom in seqnames(regions)@values){
    print(paste('Getting window-averaged methylation for ', chrom, '...', sep = ''))
    chr.index <- which(seqnames(regions) == chrom)
    tmp <- regions[chr.index, ]
    ol <- findOverlaps(tmp, AM.gr)
    srle.beta <- Rle(elementMetadata(AM.gr)$beta[subjectHits(ol)])
    srle.gamma <- Rle(elementMetadata(AM.gr)$gamma[subjectHits(ol)])
    qrle <- Rle(queryHits(ol))
    v.beta <- Views(srle.beta, successiveIRanges(runLength(qrle)))
    v.gamma <- Views(srle.gamma, successiveIRanges(runLength(qrle)))
    means$beta_w[row.index:(row.index + length(tmp) - 1)] <- viewMeans(v.beta)
    means$gamma_w[row.index:(row.index + length(tmp) - 1)] <- viewMeans(v.gamma)
    row.index <- row.index + length(tmp)
  }
  return(means)
}

#### Annotate WF_outermost.gr with CGI, n_500, beta_500 and gamma_500 ####
elementMetadata(WF_outermost.gr)$CGI <- countOverlaps(WF_outermost.gr, CGI, type = 'within')
elementMetadata(WF_outermost.gr)$n_500 <- countCpGsInWindow(WF_outermost.gr, Hsapiens, 500)
ave.meth.tmp <- averageMethylationInWindow(WF_outermost.gr, AM.gr, 500)
elementMetadata(WF_outermost.gr)$beta_500 <- ave.meth.tmp$beta_w
elementMetadata(WF_outermost.gr)$gamma_500 <- ave.meth.tmp$gamma_w
rm(ave.meth.tmp)
gc()

#### Plot aggregate log odds ratios as a function of distance, stratified by n_500 and gamma_500 ####
gamma_500.tertile <- quantile(elementMetadata(WF_outermost.gr)$gamma_500, c(1/3, 2/3))
t1 <- lor(subset(WF_outermost.gr, elementMetadata(WF_outermost.gr)$CGI == 1 & elementMetadata(WF_outermost.gr)$gamma_500 <= gamma_500.tertile[1]), strand = 'combined', aggregateByLag = TRUE)
t2 <- lor(subset(WF_outermost.gr, elementMetadata(WF_outermost.gr)$CGI == 1 & elementMetadata(WF_outermost.gr)$gamma_500 > gamma_500.tertile[1] & elementMetadata(WF_outermost.gr)$gamma_500 < gamma_500.tertile[2]), strand = 'combined', aggregateByLag = TRUE)
t3 <- lor(subset(WF_outermost.gr, elementMetadata(WF_outermost.gr)$CGI == 1 & elementMetadata(WF_outermost.gr)$gamma_500 >= gamma_500.tertile[1]), strand = 'combined', aggregateByLag = TRUE)

b1 <- lor(subset(WF_outermost.gr, elementMetadata(WF_outermost.gr)$CGI == 0 & elementMetadata(WF_outermost.gr)$gamma_500 <= gamma_500.tertile[1]), strand = 'combined', aggregateByLag = TRUE)
b2 <- lor(subset(WF_outermost.gr, elementMetadata(WF_outermost.gr)$CGI == 0 & elementMetadata(WF_outermost.gr)$gamma_500 > gamma_500.tertile[1] & elementMetadata(WF_outermost.gr)$gamma_500 < gamma_500.tertile[2]), strand = 'combined', aggregateByLag = TRUE)
b3 <- lor(subset(WF_outermost.gr, elementMetadata(WF_outermost.gr)$CGI == 0 & elementMetadata(WF_outermost.gr)$gamma_500 >= gamma_500.tertile[1]), strand = 'combined', aggregateByLag = TRUE)

lor.df <- rbind(t1, t2, t3, b1, b2, b3)
lor.df$CGI <- c(rep('CGI', nrow(t1) + nrow(t2) + nrow(t3)), rep('Non-CGI', nrow(b1) + nrow(b2) + nrow(b3)))
lor.df$gamma <- c(rep('Low methylation', nrow(t1)), rep('Intermediate methylation', nrow(t2)), rep('High methylation', nrow(t3)), rep('Low methylation', nrow(b1)), rep('Intermediate methylation', nrow(b2)), rep('High methylation', nrow(b3)))

p <- ggplot(data = lor.df, aes(x = lag, y = lor)) + geom_point() + facet_grid(CGI ~ gamma) + opts(title = paste(sample.name, ': Within-fragment aggregated comethylation')) + scale_x_continuous('Distance between CpGs (bp)') + scale_y_continuous('Log odds ratio')
ggsave(p, filename = paste('plots/', sample.name, '_lor_by_methylation_level_aggregated_by_distance_stratified.pdf', sep = ''))

#### Save workspace to file and write lor.df to file ####
save.image(paste(sample.name, '_lor_by_methylation_level_aggregated_by_distance.RData', sep = ''))
write.table(lor.df, paste(sample.name, '_lor_by_methylation_level_aggregated_by_distance.txt', sep = ''), sep = '\t', quote = FALSE, row.names = FALSE)
       
       
#### Finished ####
       

