#### Peter Hickey ####
# 28/06/2012
# Comethylation stratified by level of methylation

#### TODO ####
# Add average methylation level of CGI
# Remove dependencies on external packages where possible
# Vary the CpG-density window size to see how it affects the results
# Somehow take into account how many of the CpGs in the window actually have sufficient coverage to call a beta-value
# Incorporate line colour, as well as linetype, into quantile-plots (see plot_methylation_read_position_bias.r)
# Make the plots into a function
# Just stratify by average methylation-level (not CpG-density) - it looks like we can see the 10bp periodicity in this data!!!

#### Read in command line arguments ####
args <- commandArgs(TRUE)
# The name of the sample of interest
sample.name <- args[1]

#### Load libraries ####
library(stringr)
library(GenomicRanges)
library(plyr)
library(ggplot2)
library(doMC)
library(BSgenome)
library('BSgenome.Hsapiens.UCSC.hg18')
library(parallel)
library(reshape2)

#### Register doMC backend ####
registerDoMC(cores = 10)

#### Function definitions ####
# Output is a genomic ranges object with start = C in CpG1 and end = C in CpG2
# If chromosomes == 'all' then we use chr1, ..., chr22, chrX and chrY. We don't use chrM because it is a circular chromosome, which breaks subsetByOverlaps()
readWF <- function(sample.name, chromosomes = 'all', pairChoice = 'all'){
  # Construct a list of chromosome names that are to have lag-correlations computed.
  if(chromosomes == 'all'){
    chromosome.list <- list('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                            'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY')
    names(chromosome.list) <- chromosome.list
  } else{
    chromosome.list <- list(chromosomes)
    names(chromosome.list) <- chromosomes
  }
  WF.files <- list()
  for(chrom in names(chromosome.list)){
    if(pairChoice == 'all'){
      filename <- str_c(sample.name, '_', chrom, '.wf')
    } else if(pairChoice == 'outermost'){
      filename <- str_c(sample.name, '_', chrom, '_outermost.wf')
    }
    WF.files[[chrom]] <- read.table(filename, sep = '\t', stringsAsFactors = FALSE, header = TRUE)
  }
  WF.df <- rbind.fill(WF.files)
  WF.gr <- GRanges(seqnames = WF.df$chr, ranges = IRanges(start = WF.df$pos1, end = WF.df$pos2), counts = DataFrame(WF.df[, 4:15]))
  names(values(WF.gr)) <- str_sub(names(values(WF.gr)), 8)
  seqlengths(WF.gr)[names(chromosome.list)] <- seqlengths(Hsapiens)[names(chromosome.list)]
  rm(WF.df)
  return(WF.gr)
}

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
  if(aggregateByLag == TRUE){
    return.df <- data.frame(lag = agg.counts$lag, lor = lor, ASE = ASE)
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

#### Load WF and AM files ####
WF.gr <- readWF(sample.name, chromosomes = 'all', pairChoice = 'all')
AM <- read.table(str_c('../AM/', sample.name, '.am'), header = TRUE, stringsAsFactors = FALSE, )
AM.gr <- GRanges(seqnames = AM$chr, ranges = IRanges(start = AM$pos, end = AM$pos + 1), beta = AM$beta, gamma = AM$gamma)
rm(AM)
gc()

#### Load CpG islands (CGIs) ####
CGI <- read.table('/home/users/lab0605/hickey/Rafas-cpg-islands-hg18.txt', header = TRUE, stringsAsFactors = FALSE)
CGI <- GRanges(seqnames = CGI$chr, ranges = IRanges(start = CGI$start, end = CGI$end))
in.CGI <- countOverlaps(WF.gr, CGI, type = 'within')
rm(CGI)
gc()

#### Construct 'CpG.pairs' a GRanges object which contains all CpG-pairs, along with each pair's log-odds ratio, coverage, number of CpGs in a window of size 'w' (n_w, w = 500, 1000), the average methylation in each window (as measured by beta or gamma, beta_500, beta_1000, gamma_500, gamma_1000) and whether the CpG-pair is in a CGI.  ####
lors <- lor(WF.gr, strand = 'combined')$lor
coverage <- elementMetadata(WF.gr)$MM + elementMetadata(WF.gr)$MU + elementMetadata(WF.gr)$UM + elementMetadata(WF.gr)$UU
n_500 <- countCpGsInWindow(WF.gr, Hsapiens, 500)
#n_1000 <- countCpGsInWindow(WF.gr, Hsapiens, 1000)
average.methylation_500 <- averageMethylationInWindow(WF.gr, AM.gr, 500)
#average.methylation_1000 <- averageMethylationInWindow(WF.gr, AM.gr, 1000)
#CpG.pairs <- GRanges(seqnames = seqnames(WF.gr), ranges = ranges(WF.gr), lor = lors, coverage = coverage, n_500 = n_500, beta_500 = average.methylation_500$beta, gamma_500 = average.methylation_500$gamma, n_1000 = n_1000, beta_1000 = average.methylation_1000$beta, gamma_1000 = average.methylation_1000$gamma, CGI = in.CGI)
CpG.pairs <- GRanges(seqnames = seqnames(WF.gr), ranges = ranges(WF.gr), lor = lors, coverage = coverage, n_500 = n_500, beta_500 = average.methylation_500$beta_w, gamma_500 = average.methylation_500$gamma_w, CGI = in.CGI)
seqlengths(CpG.pairs) <- seqlengths(WF.gr)
#rm(WF.gr, AM.gr, lors, coverage, n_500, n_1000, average.methylation_500, average.methylation_1000)
rm(WF.gr, AM.gr, lors, coverage, n_500, average.methylation_500, in.CGI)
gc()

#### Create data.frame to be used in WF-plots ####
lags <- width(CpG.pairs) - 1
# Discretise n_500, beta_500 and gamma_500 
median.n <- median(elementMetadata(CpG.pairs)$n_500)
n_500 <- ifelse(elementMetadata(CpG.pairs)$n_500 >= median.n, "CpG-dense", "CpG-sparse")
beta.groups <- ifelse(elementMetadata(CpG.pairs)$beta_500 >= 0.8, "High methylation (beta)", ifelse(elementMetadata(CpG.pairs)$beta_500 <= 0.2, "Low methylation (beta)", "Intermediate methylation (beta)"))
beta_500 <- factor(beta.groups, labels = c("beta <= 0.2", expression(paste(0.2 < beta) < 0.8), "beta >= 0.8"))

gamma.tertiles <- quantile(elementMetadata(CpG.pairs)$gamma_500, c(1/3, 2/3))
gamma.groups <- ifelse(elementMetadata(CpG.pairs)$gamma_500 >= gamma.tertiles[2], "High methylation (gamma)", ifelse(elementMetadata(CpG.pairs)$gamma_500 <= gamma.tertiles[1], "Low methylation (gamma)", "Intermediate methylation (gamma)"))
gamma_500 <- factor(gamma.groups, labels = c("gamma <= gamma[(1/3)]", expression(paste(gamma[(1/3)] < gamma) < gamma[(2/3)]), "gamma >= gamma[(2/3)]"))

df.beta <- data.frame(lag = lags, n_500 = n_500, beta_500 = beta_500, coverage = elementMetadata(CpG.pairs)$coverage, CGI = elementMetadata(CpG.pairs)$CGI, lor = elementMetadata(CpG.pairs)$lor)
df.gamma <- data.frame(lag = lags, n_500 = n_500, gamma_500 = gamma_500, coverage = elementMetadata(CpG.pairs)$coverage, CGI = elementMetadata(CpG.pairs)$CGI, lor = elementMetadata(CpG.pairs)$lor)

rm(averageMethylationInWindow, beta.groups, beta_500, countCpGsInWindow, gamma.groups, gamma_500, gamma.tertiles, lags, lor, median.n, n_500, readWF)
gc()

#### Save key objects form workspace and change working directory to 'plots' directory ####
save(sample.name, CpG.pairs, df.beta, df.gamma, file = paste(sample.name, '_lor_by_methylation_level.RData', sep = ''))
setwd('plots')

#### Histogram of number of CpGs in window surrounding CpG-pairs ####
pdf(file = paste(sample.name, '_CpG-density_500_histogram.pdf', sep = ''))
hist(values(CpG.pairs)$n_500, xlab = 'Number', main = paste(sample.name, ': Number of CpGs in 500bp window around CpG-pairs', sep = ""), freq = FALSE)
dev.off()

#### Histogram of beta_500 and gamma_500 ####
pdf(file = paste(sample.name, '_beta_500_histogram.pdf', sep = ''))
hist(values(CpG.pairs)$beta_500, xlab = 'beta', main = paste(sample.name, ': Average beta-values in 500bp window around CpG-pairs', sep = ""), freq = FALSE)
dev.off()

pdf(file = paste(sample.name, '_gamma_500_histogram.pdf', sep = ''))
hist(values(CpG.pairs)$gamma_500, xlab = 'gamma', main = paste(sample.name, ': Average gamma-values in 500bp window around CpG-pairs', sep = ""), freq = FALSE)
dev.off()

#### ECDF of read-depth of CpG-pairs ####
pdf(file = paste(sample.name, '_CpG-pair_coverage_ECDF.pdf', sep = ''))
plot(ecdf(elementMetadata(CpG.pairs)$coverage), main = paste(sample.name, ': ECDF of coverage of CpG-pairs', sep = ''))
dev.off()

#### ECDF of beta_500 and gamma_500 ####
pdf(file = paste(sample.name, '_beta_500_ECDF.pdf', sep = ''))
plot(ecdf(elementMetadata(CpG.pairs)$beta_500), main = paste(sample.name, ': ECDF of beta-values in 500bp window around CpG-pairs', sep = ''))
dev.off()

pdf(file = paste(sample.name, '_gamma_500_ECDF.pdf', sep = ''))
plot(ecdf(elementMetadata(CpG.pairs)$gamma_500), main = paste(sample.name, ': ECDF of gamma-values in 500bp window around CpG-pairs', sep = ''))
dev.off()

#### All CpG-pairs: LOR as a function of intra-pair distance, stratified by CpG-density and average beta-value of region around CpG-pair ####
# Beta
df.beta.quantile.wide <- ddply(df.beta, .variables = c('lag', 'n_500', 'beta_500'), function(x){quantile(x$lor, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)}, .parallel = TRUE)
df.beta.quantile <- melt(df.beta.quantile.wide, c('lag', 'n_500', 'beta_500'), variable.name = "Quantile", value.name = "LOR")
b1 <- qplot(x = lag, y = LOR, lty = Quantile, data = df.beta.quantile, geom = "line", main = paste(sample.name, ': WF-comethylation all CpG-pairs (n = ', nrow(df.beta), ')', sep = ""), ylab = "Log odds ratio", xlab = "Distance between CpGs (bp)") + facet_grid(n_500 ~ beta_500, label = label_parsed) + scale_linetype_manual(name="Quantile",  values=c('10%' = 3, '25%' = 2, '50%' = 1, '75%' = 2, '90%' = 3))
ggsave(filename = paste(sample.name, '_lor_by_beta_level_all_CpG-pairs.pdf', sep = ''), plot = b1)

# Gamma
df.gamma.quantile.wide <- ddply(df.gamma, .variables = c('lag', 'n_500', 'gamma_500'), function(x){quantile(x$lor, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)}, .parallel = TRUE)
df.gamma.quantile <- melt(df.gamma.quantile.wide, c('lag', 'n_500', 'gamma_500'), variable.name = "Quantile", value.name = "LOR")
g1 <- qplot(x = lag, y = LOR, lty = Quantile, data = df.gamma.quantile, geom = "line", main = paste(sample.name, ': WF-comethylation all CpG-pairs (n = ', nrow(df.gamma), ')', sep = ""), ylab = "Log odds ratio", xlab = "Distance between CpGs (bp)") + facet_grid(n_500 ~ gamma_500, label = label_parsed) + scale_linetype_manual(name="Quantile",  values=c('10%' = 3, '25%' = 2, '50%' = 1, '75%' = 2, '90%' = 3))
ggsave(filename = paste(sample.name, '_lor_by_gamma_level_all_CpG-pairs.pdf', sep = ''), plot = g1)

#### All CpG-pairs with coverage >= 50% of CpG-pair-coverage: LOR as a function of intra-pair distance, stratified by CpG-density and average beta-value of region around CpG-pair ####
# Beta
min.cov.50 <- quantile(df.beta$coverage, 0.5)
df.beta.quantile.50.wide <- ddply(subset(df.beta, coverage >= min.cov.50), .variables = c('lag', 'n_500', 'beta_500'), function(x){quantile(x$lor, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)}, .parallel = TRUE)
df.beta.quantile.50 <- melt(df.beta.quantile.50.wide, c('lag', 'n_500', 'beta_500'), variable.name = "Quantile", value.name = "LOR")
b2 <- qplot(x = lag, y = LOR, lty = Quantile, data = df.beta.quantile.50, geom = "line", main = paste(sample.name, ': WF-comethylation >= 50% quantile of pair-coverage (= ', min.cov.50, 'x)', sep = ""), ylab = "log-odds ratio", xlab = "Distance between CpGs (bp)") + facet_grid(n_500 ~ beta_500, label = label_parsed) + scale_linetype_manual(name="Quantile",  values=c('10%' = 3, '25%' = 2, '50%' = 1, '75%' = 2, '90%' = 3))
ggsave(filename = paste(sample.name, '_lor_by_beta_level_CpG-pairs_with_cov_geq_Q50.pdf', sep = ''), plot = b2)
if(max(df.beta$lag) > 100){
  b2a <- b2 + scale_x_continuous(limits = c(0, 100))
  ggsave(filename = paste(sample.name, '_lor_by_beta_level_CpG-pairs_with_cov_geq_Q50_0to100.pdf', sep = ''), plot = b2a)
}

# Gamma
min.cov.50 <- quantile(df.gamma$coverage, 0.5)
df.gamma.quantile.50.wide <- ddply(subset(df.gamma, coverage >= min.cov.50), .variables = c('lag', 'n_500', 'gamma_500'), function(x){quantile(x$lor, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)}, .parallel = TRUE)
df.gamma.quantile.50 <- melt(df.gamma.quantile.50.wide, c('lag', 'n_500', 'gamma_500'), variable.name = "Quantile", value.name = "LOR")
g2 <- qplot(x = lag, y = LOR, lty = Quantile, data = df.gamma.quantile.50, geom = "line", main = paste(sample.name, ': WF-comethylation >= 50% quantile of pair-coverage (= ', min.cov.50, 'x)', sep = ""), ylab = "log-odds ratio", xlab = "Distance between CpGs (bp)") + facet_grid(n_500 ~ gamma_500, label = label_parsed) + scale_linetype_manual(name="Quantile",  values=c('10%' = 3, '25%' = 2, '50%' = 1, '75%' = 2, '90%' = 3))
ggsave(filename = paste(sample.name, '_lor_by_gamma_level_CpG-pairs_with_cov_geq_Q50.pdf', sep = ''), plot = g2)
if(max(df.gamma$lag) > 100){
  g2a <- g2 + scale_x_continuous(limits = c(0, 100))
  ggsave(filename = paste(sample.name, '_lor_by_gamma_level_CpG-pairs_with_cov_geq_Q50_0to100.pdf', sep = ''), plot = g2a)
}

#### All CpG-pairs with coverage >= 75% of CpG-pair-coverage: LOR as a function of intra-pair distance, stratified by CpG-density and average beta-value of region around CpG-pair ####
# Beta
min.cov.75 <- quantile(df.beta$coverage, 0.75)
df.beta.quantile.75.wide <- ddply(subset(df.beta, coverage >= min.cov.75), .variables = c('lag', 'n_500', 'beta_500'), function(x){quantile(x$lor, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)}, .parallel = TRUE)
df.beta.quantile.75 <- melt(df.beta.quantile.75.wide, c('lag', 'n_500', 'beta_500'), variable.name = "Quantile", value.name = "LOR")
b3 <- qplot(x = lag, y = LOR, lty = Quantile, data = df.beta.quantile.75, geom = "line", main = paste(sample.name, ': WF-comethylation >= 75% quantile of pair-coverage (= ', min.cov.75, 'x)', sep = ""), ylab = "log-odds ratio", xlab = "Distance between CpGs (bp)") + facet_grid(n_500 ~ beta_500, label = label_parsed) + scale_linetype_manual(name="Quantile",  values=c('10%' = 3, '25%' = 2, '50%' = 1, '75%' = 2, '90%' = 3))
ggsave(filename = paste(sample.name, '_lor_by_beta_level_CpG-pairs_with_cov_geq_Q75.pdf', sep = ''), plot = b3)
if(max(df.beta$lag) > 100){
  b3a <- b3 + scale_x_continuous(limits = c(0, 100))
  ggsave(filename = paste(sample.name, '_lor_by_beta_level_CpG-pairs_with_cov_geq_Q75_0to100.pdf', sep = ''), plot = b3a)
}

# Gamma
min.cov.75 <- quantile(df.gamma$coverage, 0.75)
df.gamma.quantile.75.wide <- ddply(subset(df.gamma, coverage >= min.cov.75), .variables = c('lag', 'n_500', 'gamma_500'), function(x){quantile(x$lor, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)}, .parallel = TRUE)
df.gamma.quantile.75 <- melt(df.gamma.quantile.75.wide, c('lag', 'n_500', 'gamma_500'), variable.name = "Quantile", value.name = "LOR")
g3 <- qplot(x = lag, y = LOR, lty = Quantile, data = df.gamma.quantile.75, geom = "line", main = paste(sample.name, ': WF-comethylation >= 75% quantile of pair-coverage (= ', min.cov.75, 'x)', sep = ""), ylab = "log-odds ratio", xlab = "Distance between CpGs (bp)") + facet_grid(n_500 ~ gamma_500, label = label_parsed) + scale_linetype_manual(name="Quantile",  values=c('10%' = 3, '25%' = 2, '50%' = 1, '75%' = 2, '90%' = 3))
ggsave(filename = paste(sample.name, '_lor_by_gamma_level_CpG-pairs_with_cov_geq_Q75.pdf', sep = ''), plot = g3)
if(max(df.gamma$lag) > 100){
  g3a <- g3 + scale_x_continuous(limits = c(0, 100))
  ggsave(filename = paste(sample.name, '_lor_by_gamma_level_CpG-pairs_with_cov_geq_Q75_0to100.pdf', sep = ''), plot = g3a)
}

#### All CpG-pairs with coverage >= 90% of CpG-pair-coverage: LOR as a function of intra-pair distance, stratified by CpG-density and average beta-value of region around CpG-pair ####
# Beta
min.cov.90 <- quantile(df.beta$coverage, 0.90)
df.beta.quantile.90.wide <- ddply(subset(df.beta, coverage >= min.cov.90), .variables = c('lag', 'n_500', 'beta_500'), function(x){quantile(x$lor, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)}, .parallel = TRUE)
df.beta.quantile.90 <- melt(df.beta.quantile.90.wide, c('lag', 'n_500', 'beta_500'), variable.name = "Quantile", value.name = "LOR")
b4 <- qplot(x = lag, y = LOR, lty = Quantile, data = df.beta.quantile.90, geom = "line", main = paste(sample.name, ': WF-comethylation >= 90% quantile of pair-coverage (= ', min.cov.90, 'x)', sep = ""), ylab = "log-odds ratio", xlab = "Distance between CpGs (bp)") + facet_grid(n_500 ~ beta_500, label = label_parsed) + scale_linetype_manual(name="Quantile",  values=c('10%' = 3, '25%' = 2, '50%' = 1, '75%' = 2, '90%' = 3))
ggsave(filename = paste(sample.name, '_lor_by_beta_level_CpG-pairs_with_cov_geq_Q90.pdf', sep = ''), plot = b4)
if(max(df.beta$lag) > 100){
  b4a <- b4 + scale_x_continuous(limits = c(0, 100))
  ggsave(filename = paste(sample.name, '_lor_by_beta_level_CpG-pairs_with_cov_geq_Q90_0to100.pdf', sep = ''), plot = b4a)
}

# Gamma
min.cov.90 <- quantile(df.gamma$coverage, 0.90)
df.gamma.quantile.90.wide <- ddply(subset(df.gamma, coverage >= min.cov.90), .variables = c('lag', 'n_500', 'gamma_500'), function(x){quantile(x$lor, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)}, .parallel = TRUE)
df.gamma.quantile.90 <- melt(df.gamma.quantile.90.wide, c('lag', 'n_500', 'gamma_500'), variable.name = "Quantile", value.name = "LOR")
g4 <- qplot(x = lag, y = LOR, lty = Quantile, data = df.gamma.quantile.90, geom = "line", main = paste(sample.name, ': WF-comethylation >= 90% quantile of pair-coverage (= ', min.cov.90, 'x)', sep = ""), ylab = "log-odds ratio", xlab = "Distance between CpGs (bp)") + facet_grid(n_500 ~ gamma_500, label = label_parsed) + scale_linetype_manual(name="Quantile",  values=c('10%' = 3, '25%' = 2, '50%' = 1, '75%' = 2, '90%' = 3))
ggsave(filename = paste(sample.name, '_lor_by_gamma_level_CpG-pairs_with_cov_geq_Q90.pdf', sep = ''), plot = g4)
if(max(df.gamma$lag) > 100){
  g4a <- g4 + scale_x_continuous(limits = c(0, 100))
  ggsave(filename = paste(sample.name, '_lor_by_gamma_level_CpG-pairs_with_cov_geq_Q90_0to100.pdf', sep = ''), plot = g4a)
}

#### All CpG-pairs: LOR as a function of intra-pair distance, stratified by CGI-status and average beta-value of region around CpG-pair ####
# Beta
df.beta.quantile.wide <- ddply(df.beta, .variables = c('lag', 'CGI', 'beta_500'), function(x){quantile(x$lor, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)}, .parallel = TRUE)
df.beta.quantile <- melt(df.beta.quantile.wide, c('lag', 'CGI', 'beta_500'), variable.name = "Quantile", value.name = "LOR")
b1_CGI <- qplot(x = lag, y = LOR, lty = Quantile, data = df.beta.quantile, geom = "line", main = paste(sample.name, ': WF-comethylation all CpG-pairs (n = ', nrow(df.beta), ')', sep = ""), ylab = "Log odds ratio", xlab = "Distance between CpGs (bp)") + facet_grid(CGI ~ beta_500, label = label_parsed) + scale_linetype_manual(name="Quantile",  values=c('10%' = 3, '25%' = 2, '50%' = 1, '75%' = 2, '90%' = 3))
ggsave(filename = paste(sample.name, '_lor_by_beta_level_all_CpG-pairs_CGI_stratification.pdf', sep = ''), plot = b1_CGI)

# Gamma
df.gamma.quantile.wide <- ddply(df.gamma, .variables = c('lag', 'CGI', 'gamma_500'), function(x){quantile(x$lor, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)}, .parallel = TRUE)
df.gamma.quantile <- melt(df.gamma.quantile.wide, c('lag', 'CGI', 'gamma_500'), variable.name = "Quantile", value.name = "LOR")
g1_CGI <- qplot(x = lag, y = LOR, lty = Quantile, data = df.gamma.quantile, geom = "line", main = paste(sample.name, ': WF-comethylation all CpG-pairs (n = ', nrow(df.gamma), ')', sep = ""), ylab = "Log odds ratio", xlab = "Distance between CpGs (bp)") + facet_grid(CGI ~ gamma_500, label = label_parsed) + scale_linetype_manual(name="Quantile",  values=c('10%' = 3, '25%' = 2, '50%' = 1, '75%' = 2, '90%' = 3))
ggsave(filename = paste(sample.name, '_lor_by_gamma_level_all_CpG-pairs_CGI_stratification.pdf', sep = ''), plot = g1_CGI)

#### All CpG-pairs with coverage >= 50% of CpG-pair-coverage: LOR as a function of intra-pair distance, stratified by CpG-density and average beta-value of region around CpG-pair ####
# Beta
min.cov.50 <- quantile(df.beta$coverage, 0.5)
df.beta.quantile.50.wide <- ddply(subset(df.beta, coverage >= min.cov.50), .variables = c('lag', 'CGI', 'beta_500'), function(x){quantile(x$lor, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)}, .parallel = TRUE)
df.beta.quantile.50 <- melt(df.beta.quantile.50.wide, c('lag', 'CGI', 'beta_500'), variable.name = "Quantile", value.name = "LOR")
b2_CGI <- qplot(x = lag, y = LOR, lty = Quantile, data = df.beta.quantile.50, geom = "line", main = paste(sample.name, ': WF-comethylation >= 50% quantile of pair-coverage (= ', min.cov.50, 'x)', sep = ""), ylab = "log-odds ratio", xlab = "Distance between CpGs (bp)") + facet_grid(CGI ~ beta_500, label = label_parsed) + scale_linetype_manual(name="Quantile",  values=c('10%' = 3, '25%' = 2, '50%' = 1, '75%' = 2, '90%' = 3))
ggsave(filename = paste(sample.name, '_lor_by_beta_level_CpG-pairs_with_cov_geq_Q50_CGI_stratification.pdf', sep = ''), plot = b2_CGI)
if(max(df.beta$lag) > 100){
  b2a_CGI <- b2_CGI + scale_x_continuous(limits = c(0, 100))
  ggsave(filename = paste(sample.name, '_lor_by_beta_level_CpG-pairs_with_cov_geq_Q50_0to100_CGI_stratification.pdf', sep = ''), plot = b2a_CGI)
}

# Gamma
min.cov.50 <- quantile(df.gamma$coverage, 0.5)
df.gamma.quantile.50.wide <- ddply(subset(df.gamma, coverage >= min.cov.50), .variables = c('lag', 'CGI', 'gamma_500'), function(x){quantile(x$lor, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)}, .parallel = TRUE)
df.gamma.quantile.50 <- melt(df.gamma.quantile.50.wide, c('lag', 'CGI', 'gamma_500'), variable.name = "Quantile", value.name = "LOR")
g2_CGI <- qplot(x = lag, y = LOR, lty = Quantile, data = df.gamma.quantile.50, geom = "line", main = paste(sample.name, ': WF-comethylation >= 50% quantile of pair-coverage (= ', min.cov.50, 'x)', sep = ""), ylab = "log-odds ratio", xlab = "Distance between CpGs (bp)") + facet_grid(CGI ~ gamma_500, label = label_parsed) + scale_linetype_manual(name="Quantile",  values=c('10%' = 3, '25%' = 2, '50%' = 1, '75%' = 2, '90%' = 3))
ggsave(filename = paste(sample.name, '_lor_by_gamma_level_CpG-pairs_with_cov_geq_Q50_CGI_stratification.pdf', sep = ''), plot = g2_CGI)
if(max(df.gamma$lag) > 100){
  g2a_CGI <- g2_CGI + scale_x_continuous(limits = c(0, 100))
  ggsave(filename = paste(sample.name, '_lor_by_gamma_level_CpG-pairs_with_cov_geq_Q50_0to100_CGI_stratification.pdf', sep = ''), plot = g2a_CGI)
}

#### All CpG-pairs with coverage >= 75% of CpG-pair-coverage: LOR as a function of intra-pair distance, stratified by CpG-density and average beta-value of region around CpG-pair ####
# Beta
min.cov.75 <- quantile(df.beta$coverage, 0.75)
df.beta.quantile.75.wide <- ddply(subset(df.beta, coverage >= min.cov.75), .variables = c('lag', 'CGI', 'beta_500'), function(x){quantile(x$lor, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)}, .parallel = TRUE)
df.beta.quantile.75 <- melt(df.beta.quantile.75.wide, c('lag', 'CGI', 'beta_500'), variable.name = "Quantile", value.name = "LOR")
b3_CGI <- qplot(x = lag, y = LOR, lty = Quantile, data = df.beta.quantile.75, geom = "line", main = paste(sample.name, ': WF-comethylation >= 75% quantile of pair-coverage (= ', min.cov.75, 'x)', sep = ""), ylab = "log-odds ratio", xlab = "Distance between CpGs (bp)") + facet_grid(CGI ~ beta_500, label = label_parsed) + scale_linetype_manual(name="Quantile",  values=c('10%' = 3, '25%' = 2, '50%' = 1, '75%' = 2, '90%' = 3))
ggsave(filename = paste(sample.name, '_lor_by_beta_level_CpG-pairs_with_cov_geq_Q75_CGI_stratification.pdf', sep = ''), plot = b3_CGI)
if(max(df.beta$lag) > 100){
  b3a_CGI <- b3_CGI + scale_x_continuous(limits = c(0, 100))
  ggsave(filename = paste(sample.name, '_lor_by_beta_level_CpG-pairs_with_cov_geq_Q75_0to100_CGI_stratification.pdf', sep = ''), plot = b3a_CGI)
}

# Gamma
min.cov.75 <- quantile(df.gamma$coverage, 0.75)
df.gamma.quantile.75.wide <- ddply(subset(df.gamma, coverage >= min.cov.75), .variables = c('lag', 'CGI', 'gamma_500'), function(x){quantile(x$lor, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)}, .parallel = TRUE)
df.gamma.quantile.75 <- melt(df.gamma.quantile.75.wide, c('lag', 'CGI', 'gamma_500'), variable.name = "Quantile", value.name = "LOR")
g3_CGI <- qplot(x = lag, y = LOR, lty = Quantile, data = df.gamma.quantile.75, geom = "line", main = paste(sample.name, ': WF-comethylation >= 75% quantile of pair-coverage (= ', min.cov.75, 'x)', sep = ""), ylab = "log-odds ratio", xlab = "Distance between CpGs (bp)") + facet_grid(CGI ~ gamma_500, label = label_parsed) + scale_linetype_manual(name="Quantile",  values=c('10%' = 3, '25%' = 2, '50%' = 1, '75%' = 2, '90%' = 3))
ggsave(filename = paste(sample.name, '_lor_by_gamma_level_CpG-pairs_with_cov_geq_Q75_CGI_stratification.pdf', sep = ''), plot = g3_CGI)
if(max(df.gamma$lag) > 100){
  g3a_CGI <- g3_CGI + scale_x_continuous(limits = c(0, 100))
  ggsave(filename = paste(sample.name, '_lor_by_gamma_level_CpG-pairs_with_cov_geq_Q75_0to100_CGI_stratification.pdf', sep = ''), plot = g3a_CGI)
}

#### All CpG-pairs with coverage >= 90% of CpG-pair-coverage: LOR as a function of intra-pair distance, stratified by CpG-density and average beta-value of region around CpG-pair ####
# Beta
min.cov.90 <- quantile(df.beta$coverage, 0.90)
df.beta.quantile.90.wide <- ddply(subset(df.beta, coverage >= min.cov.90), .variables = c('lag', 'CGI', 'beta_500'), function(x){quantile(x$lor, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)}, .parallel = TRUE)
df.beta.quantile.90 <- melt(df.beta.quantile.90.wide, c('lag', 'CGI', 'beta_500'), variable.name = "Quantile", value.name = "LOR")
b4_CGI <- qplot(x = lag, y = LOR, lty = Quantile, data = df.beta.quantile.90, geom = "line", main = paste(sample.name, ': WF-comethylation >= 90% quantile of pair-coverage (= ', min.cov.90, 'x)', sep = ""), ylab = "log-odds ratio", xlab = "Distance between CpGs (bp)") + facet_grid(CGI ~ beta_500, label = label_parsed) + scale_linetype_manual(name="Quantile",  values=c('10%' = 3, '25%' = 2, '50%' = 1, '75%' = 2, '90%' = 3))
ggsave(filename = paste(sample.name, '_lor_by_beta_level_CpG-pairs_with_cov_geq_Q90_CGI_stratification.pdf', sep = ''), plot = b4_CGI)
if(max(df.beta$lag) > 100){
  b4a_CGI <- b4_CGI + scale_x_continuous(limits = c(0, 100))
  ggsave(filename = paste(sample.name, '_lor_by_beta_level_CpG-pairs_with_cov_geq_Q90_0to100_CGI_stratification.pdf', sep = ''), plot = b4a_CGI)
}

# Gamma
min.cov.90 <- quantile(df.gamma$coverage, 0.90)
df.gamma.quantile.90.wide <- ddply(subset(df.gamma, coverage >= min.cov.90), .variables = c('lag', 'CGI', 'gamma_500'), function(x){quantile(x$lor, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)}, .parallel = TRUE)
df.gamma.quantile.90 <- melt(df.gamma.quantile.90.wide, c('lag', 'CGI', 'gamma_500'), variable.name = "Quantile", value.name = "LOR")
g4_CGI <- qplot(x = lag, y = LOR, lty = Quantile, data = df.gamma.quantile.90, geom = "line", main = paste(sample.name, ': WF-comethylation >= 90% quantile of pair-coverage (= ', min.cov.90, 'x)', sep = ""), ylab = "log-odds ratio", xlab = "Distance between CpGs (bp)") + facet_grid(CGI ~ gamma_500, label = label_parsed) + scale_linetype_manual(name="Quantile",  values=c('10%' = 3, '25%' = 2, '50%' = 1, '75%' = 2, '90%' = 3))
ggsave(filename = paste(sample.name, '_lor_by_gamma_level_CpG-pairs_with_cov_geq_Q90_CGI_stratification.pdf', sep = ''), plot = g4_CGI)
if(max(df.gamma$lag) > 100){
  g4a_CGI <- g4_CGI + scale_x_continuous(limits = c(0, 100))
  ggsave(filename = paste(sample.name, '_lor_by_gamma_level_CpG-pairs_with_cov_geq_Q90_0to100_CGI_stratification.pdf', sep = ''), plot = g4a_CGI)
}

#### Finished ####

