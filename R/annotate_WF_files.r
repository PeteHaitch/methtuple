# Peter Hickey
# 10/09/2012
# Annotate WF files with coverage, log odds ratios, NIC etc.
# 

#### TODOs ####
# Implement a "WF" class
# Make all functions act on "WF" objects
# Remove genomic-features if not required

#### Read in command line arguments ####
args <- commandArgs(TRUE)
# The name of the sample of interest
sample.name <- args[1]
# Path to plot_comethylation_v2_functions.r
path.to.functions <- args[2]

#### Load libraries ####
library(rtracklayer) # I use rtracklayer::import to get the reference CpGs from a bed file

#### Import hg18 reference CpGs ####
cpgs <- import('~/CpGs_hg18.bed.gz', asRangedData = FALSE)

#### Source functions ####
source(path.to.functions)

#### Read in the data ####
wf.all.gr <- read.wf(sample.name, chromosomes = 'all', pair.choice = 'all')
wf.outermost.gr <- read.wf(sample.name, chromosomes = 'all', pair.choice = 'outermost')

#### Compute log-odds ratios for each CpG-pair ####
tmp <- lor(values(wf.all.gr), correct = TRUE)
elementMetadata(wf.all.gr)$lor <- tmp$lor
elementMetadata(wf.all.gr)$ase <- tmp$ase
tmp <- lor(values(wf.outermost.gr), correct = TRUE)
elementMetadata(wf.outermost.gr)$lor <- tmp$lor
elementMetadata(wf.outermost.gr)$ase <- tmp$ase

#### Compute coverage of each CpG-pair ####
elementMetadata(wf.all.gr)$cov <- elementMetadata(wf.all.gr)$MM + elementMetadata(wf.all.gr)$MU + elementMetadata(wf.all.gr)$UM + elementMetadata(wf.all.gr)$UU
elementMetadata(wf.outermost.gr)$cov <- elementMetadata(wf.outermost.gr)$MM + elementMetadata(wf.outermost.gr)$MU + elementMetadata(wf.outermost.gr)$UM + elementMetadata(wf.outermost.gr)$UU

#### Compute number of intervening CpGs (NIC) for each CpG-pair ####
elementMetadata(wf.all.gr)$NIC <- countOverlaps(wf.all.gr, cpgs) - 2 # All CpG-pairs overlap 2 CpGs by definition; therefore need to subtract 2 to get NIC
elementMetadata(wf.outermost.gr)$NIC <- countOverlaps(wf.outermost.gr, cpgs) - 2 # All CpG-pairs overlap 2 CpGs by definition; therefore need to subtract 2 to get NIC

#### Create WF objects (both all and outermost) where all CpG-pairs have NIC==0 ####
wf.zero.nic.all.gr <- wf.all.gr[values(wf.all.gr)$NIC == 0, ]
wf.zero.nic.outermost.gr <- wf.outermost.gr[values(wf.outermost.gr)$NIC == 0, ]

#### Load AM files ####
am <- read.table(paste0('../AM/', sample.name, '.am'), header = TRUE, stringsAsFactors = FALSE, )
am.gr <- GRanges(seqnames = am$chr, ranges = IRanges(start = am$pos, end = am$pos + 1), beta = am$beta, gamma = am$gamma)
rm(am)
gc()

#### Annotate each CpG-pair with the average beta- and gamma-values in a 500 bp window ####
tmp <- averageMethylationInWindow(wf.all.gr, am.gr, 500)
elementMetadata(wf.all.gr)$beta.500 <- tmp$beta_w
elementMetadata(wf.all.gr)$gamma.500 <- tmp$gamma_w
rm(tmp)
tmp <- averageMethylationInWindow(wf.outermost.gr, am.gr, 500)
elementMetadata(wf.outermost.gr)$beta.500 <- tmp$beta_w
elementMetadata(wf.outermost.gr)$gamma.500 <- tmp$gamma_w
rm(tmp)
tmp <- averageMethylationInWindow(wf.zero.nic.all.gr, am.gr, 500)
elementMetadata(wf.zero.nic.all.gr)$beta.500 <- tmp$beta_w
elementMetadata(wf.zero.nic.all.gr)$gamma.500 <- tmp$gamma_w
rm(tmp)
tmp <- averageMethylationInWindow(wf.zero.nic.outermost.gr, am.gr, 500)
elementMetadata(wf.zero.nic.outermost.gr)$beta.500 <- tmp$beta_w
elementMetadata(wf.zero.nic.outermost.gr)$gamma.500 <- tmp$gamma_w
rm(tmp)

#### Save complete.wf.all.gr, complete.wf.outermost.gr, complete.wf.zero.nic.all.gr and complete.wf.zero.nic.outermost.gr as an RData object ####
complete.wf.all.gr <- wf.all.gr
complete.wf.outermost.gr <- wf.outermost.gr
complete.wf.zero.nic.all.gr <- wf.zero.nic.all.gr
complete.wf.zero.nic.outermost.gr <- wf.zero.nic.outermost.gr
save(sample.name, complete.wf.all.gr, complete.wf.outermost.gr, complete.wf.zero.nic.all.gr, complete.wf.zero.nic.outermost.gr, file = paste0(sample.name, '_complete_annotated_WF_objects.RData'))

#### Save only pair-specific variables in order to save significant space. WARNING: Aggregate log odds ratios can't be constructed from these objects. ####
values(wf.all.gr) <- values(complete.wf.all.gr)[c('lor', 'ase', 'cov', 'NIC')]
values(wf.outermost.gr) <- values(complete.wf.outermost.gr)[c('lor', 'ase', 'cov', 'NIC')]
values(wf.zero.nic.all.gr) <- values(wf.zero.nic.all.gr)[c('lor', 'ase', 'cov', 'NIC')]
values(wf.zero.nic.outermost.gr) <- values(wf.zero.nic.outermost.gr)[c('lor', 'ase', 'cov', 'NIC')]
save(sample.name, wf.all.gr, wf.outermost.gr, wf.zero.nic.all.gr, wf.zero.nic.outermost.gr, file = paste0(sample.name, '_annotated_WF_objects.RData'))

#### Finished ####


