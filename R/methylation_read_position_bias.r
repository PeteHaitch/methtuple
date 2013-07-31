#### DESCRIPTION ####
# Functions to analyse methylation read-position bias in a Bismark SAM/BAM file
# Peter Hickey (peter.hickey@gmail.com)
# 03/07/2013

#### TODO ####
# Enable separate plots for mate1 and mate2 of paired-end reads

#### LIMITATIONS ####
# Assumes no 5' trimming of reads
# Does not allow indels in reads
# directional libraries only

#### Load packages ####
library(Rsamtools)
library(stringr)

#### TESTING: parameter values ####
setwd("/stornext/User/data/allstaff/h/hickey/Seisenberger_data/J1_1/Downloaded_mapped_reads/BAM")
seq_name <- 19
seq_length <- 61431566
bam_file <- 'tmp_19.bam'
simpleCigar <- FALSE

#### Function to extract XM-tag information and summarise ####
# if simpleCigar = TRUE then only uses reads without indels
# currently does not filter any reads based on FLAG
# TODO: Add option to use a subsample of the reads (of a given size/percentage).
parse_XM <- function(seq_name, seq_length, bam_file, simpleCigar = TRUE){
  ## Construct ScanBamParam
  if (simpleCigar){
    param <- ScanBamParam(what = c('qwidth'), tag = c('XM', 'XR', 'XG'), which = GRanges(seq_name, IRanges(1, seq_length)))
  } else{
    param <- ScanBamParam(what = c('qwidth', 'cigar'), tag = c('XM', 'XR', 'XG'), which = GRanges(seq_name, IRanges(1, seq_length)))
  }
  ## Read in BAM using scanBam
  bam <- scanBam(bamFile, param = param)[[1]]
  
  ## TODO: Remove [1:100000] condition. Either use all reads or a subsample.
  ## Summarise position-specific methylation. 
  Z_count <- lapply(str_locate_all(bam$tag$XM[1:100000], pattern = 'Z'), 
                    function(x){stopifnot(x[, 'start'] == x[, 'end'])
                                x[, 'start']})
  z_count <- lapply(str_locate_all(bam$tag$XM[1:100000], pattern = 'z'), 
                    function(x){stopifnot(x[, 'start'] == x[, 'end'])
                                x[, 'start']})
  X_count <- lapply(str_locate_all(bam$tag$XM[1:100000], pattern = 'X'), 
                    function(x){stopifnot(x[, 'start'] == x[, 'end'])
                                x[, 'start']})
  x_count <- lapply(str_locate_all(bam$tag$XM[1:100000], pattern = 'x'), 
                    function(x){stopifnot(x[, 'start'] == x[, 'end'])
                                x[, 'start']})
  H_count <- lapply(str_locate_all(bam$tag$XM[1:100000], pattern = 'H'), 
                    function(x){stopifnot(x[, 'start'] == x[, 'end'])
                                x[, 'start']})
  h_count <- lapply(str_locate_all(bam$tag$XM[1:100000], pattern = 'h'), 
                    function(x){stopifnot(x[, 'start'] == x[, 'end'])
                                x[, 'start']})
  U_count <- lapply(str_locate_all(bam$tag$XM[1:100000], pattern = 'U'), 
                    function(x){stopifnot(x[, 'start'] == x[, 'end'])
                                x[, 'start']})
  u_count <- lapply(str_locate_all(bam$tag$XM[1:100000], pattern = 'u'), 
                    function(x){stopifnot(x[, 'start'] == x[, 'end'])
                                x[, 'start']})
  ## Stratify *_count by qwidth, XR and XR, and then tabulate counts to produce *_table. *_table is an 3-dimensional array of tables. The dimensions of the array are as follows: (1) qwidth; (2) XR; (3) XG; thus, for example, Z_table['90', 'CT', 'CT'] is the table of methylated CpGs by read position in reads with qwidth = 90, XR = CT and XG = CT.
  Z_table <- tapply(X = Z_count, INDEX = list(bam$qwidth[1:100000], bam$tag$XR[1:100000], bam$tag$XG[1:100000]), FUN = function(x){table(c(x, recursive = TRUE))}, simplify = TRUE)
  z_table <- tapply(X = z_count, INDEX = list(bam$qwidth[1:100000], bam$tag$XR[1:100000], bam$tag$XG[1:100000]), FUN = function(x){table(c(x, recursive = TRUE))}, simplify = TRUE)
  X_table <- tapply(X = X_count, INDEX = list(bam$qwidth[1:100000], bam$tag$XR[1:100000], bam$tag$XG[1:100000]), FUN = function(x){table(c(x, recursive = TRUE))}, simplify = TRUE)
  x_table <- tapply(X = x_count, INDEX = list(bam$qwidth[1:100000], bam$tag$XR[1:100000], bam$tag$XG[1:100000]), FUN = function(x){table(c(x, recursive = TRUE))}, simplify = TRUE)
  H_table <- tapply(X = H_count, INDEX = list(bam$qwidth[1:100000], bam$tag$XR[1:100000], bam$tag$XG[1:100000]), FUN = function(x){table(c(x, recursive = TRUE))}, simplify = TRUE)
  h_table <- tapply(X = h_count, INDEX = list(bam$qwidth[1:100000], bam$tag$XR[1:100000], bam$tag$XG[1:100000]), FUN = function(x){table(c(x, recursive = TRUE))}, simplify = TRUE)
  U_table <- tapply(X = U_count, INDEX = list(bam$qwidth[1:100000], bam$tag$XR[1:100000], bam$tag$XG[1:100000]), FUN = function(x){table(c(x, recursive = TRUE))}, simplify = TRUE)
  u_table <- tapply(X = u_count, INDEX = list(bam$qwidth[1:100000], bam$tag$XR[1:100000], bam$tag$XG[1:100000]), FUN = function(x){table(c(x, recursive = TRUE))}, simplify = TRUE)

  ## TODO: Appropriately combine counts across XR- and XG-tags to reflect their common orientation
}
