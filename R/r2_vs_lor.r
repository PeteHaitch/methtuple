# Peter Hickey
# 02/07/2012
# A quick comparison of r^2 vs. log odds ratios for measuring within-fragment comethylation

#### Load libraries ####
library(stringr)
library(GenomicRanges)
library(plyr)
library(ggplot2)
library(doMC)
library(BSgenome)
library('BSgenome.Hsapiens.UCSC.hg18')

#### Register the doMC backend ####
registerDoMC(cores = n.cores)

#### Load data ####
load('IMR90_r1_comethylation.RData')

# Plot aggregated r2 computed from WF_outermost.gr
plotAggregatedr2ByLag <- function(x, title, error.bars = TRUE, ...){
  ggplot(data = x, aes(x = lag, y = r2)) + geom_point() + opts(title = title) + scale_x_continuous('Distance between CpGs (bp)') + scale_y_continuous(limits = c(0, 1), expression(r^2)) 
}

# Plot aggregated LORs computed from WF_outermost.gr
plotAggregatedLorByLag <- function(x, title, error.bars = TRUE, ...){
  if(error.bars == TRUE){
    error.bars <- aes(ymax = lor + 2*ASE, ymin = lor - 2*ASE)
    ggplot(data = x, aes(x = lag, y = lor)) + geom_point() + opts(title = title) + scale_x_continuous('Distance between CpGs (bp)') + scale_y_continuous(limits = c(-2, 7), 'log-odds ratio') + geom_errorbar(error.bars)
  } else{
    ggplot(data = x, aes(x = lag, y = lor)) + geom_point() + opts(title = title) + scale_x_continuous('Distance between CpGs (bp)') + scale_y_continuous(limits = c(-2, 7), 'log-odds ratio') 
  }
}

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

r2 <- function(x, strand = 'combined', aggregateByLag = FALSE, aggregateByNIC = FALSE){
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
  r2 <- ifelse(MM == 0.5 & MU == 0.5 & UM == 0.5 & UU == 0.5, NA, (MM * UU - MU * UM)^2/((MM + MU) * (UM + UU) * (MU + UU) * (MM + UM)))
  if(aggregateByLag == TRUE){
    return.df <- data.frame(lag = agg.counts$lag, r2 = r2)
  } else{
    return.df <- data.frame(r2)
  }
  return(return.df)
}

a <- r2(WF_outermost.gr, aggregateByLag = TRUE)
b <- lor(WF_outermost.gr, aggregateByLag = TRUE)

plotAggregatedr2ByLag(a, expression(paste('IMR90_r1: Within-fragment comethylation using ', r^2, sep = '')), error.bars = FALSE)
ggsave('~/IMR90_r1_rsquared.pdf')
plotAggregatedLorByLag(b, 'IMR90_r1: Within-fragment comethylation using log odds ratios', error.bars = FALSE)
ggsave('~/IMR90_r1_LOR.pdf')
