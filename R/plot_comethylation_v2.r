# Peter Hickey
# 10/09/2012
# Plots of within-fragment comethylation
# This is version 2 and is a substantial re-write of the existing code.
# Functions are stored separated in plot_comethylation_v2_functions.r

#### TODOs ####
# Implement a "WF" class
# Make all functions act on "WF" objects

#### Read in command line arguments ####
args <- commandArgs(TRUE)
# The name of the sample of interest
sample.name <- args[1]
# Path to genomic-features.RData object
path.to.gf <- args[2]
# Path to plot_comethylation_v2_functions.r
path.to.functions <- args[3]

#### Load libraries ####
library(GenomicRanges)
library(ggplot2)
library(rtracklayer)
library(reshape2)
library(plyr) # I use plyr::rbind.fill in read.wf
library(BSgenome)
library('BSgenome.Hsapiens.UCSC.hg18') # Used to add seqlengths to the WF objects

#### Load genomic-features.RData, created by running create_genomic_features ####
load(file = path.to.gf)

#### Import hg18 reference CpGs ####
cpgs <- import('~/CpGs_hg18.bed.gz', asRangedData = FALSE)

#### Source functions ####
source(path.to.functions)

#### Read in the data ####
WF.all.gr <- read.wf(sample.name, chromosomes = 'all', pair.choice = 'all')
WF.outermost.gr <- read.wf(sample.name, chromosomes = 'all', pair.choice = 'outermost')

#### Compute log-odds ratios for each CpG-pair ####
tmp <- lor(values(WF.gr), correct = TRUE)
elementMetadata(WF.gr)$lor <- tmp$lor
elementMetadata(WF.gr)$ase <- tmp$ase
tmp <- lor(values(WF.outermost.gr), correct = TRUE)
elementMetadata(WF.outermost.gr)$lor <- tmp$lor
elementMetadata(WF.outermost.gr)$ase <- tmp$ase

#### Compute coverage of each CpG-pair ####
elementMetadata(WF.gr)$cov <- elementMetadata(WF.gr)$MM + elementMetadata(WF.gr)$MU + elementMetadata(WF.gr)$UM + elementMetadata(WF.gr)$UU
elementMetadata(WF.outermost.gr)$cov <- elementMetadata(WF.outermost.gr)$MM + elementMetadata(WF.outermost.gr)$MU + elementMetadata(WF.outermost.gr)$UM + elementMetadata(WF.outermost.gr)$UU

#### Compute number of intervening CpGs (NIC) for each CpG-pair
elementMetadata(WF.gr)$NIC <- countOverlaps(WF.gr, cpgs) - 2 # All CpG-pairs overlap 2 CpGs by definition; therefore need to subtract 2 to get NIC
elementMetadata(WF.outermost.gr)$NIC <- countOverlaps(WF.gr, cpgs) - 2 # All CpG-pairs overlap 2 CpGs by definition; therefore need to subtract 2 to get NIC