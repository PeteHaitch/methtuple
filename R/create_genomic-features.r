# Peter Hickey
# 01/06/2012
# Create genomic-features which will be used to stratify plots of within-fragment comethylation

#### TODOs ####
# TODO: coding regions (cds) of knownGene
# TODO: Find out whether knownGene is the best database to use

#### Load libraries ####
library(stringr)
library(GenomicRanges)
library(plyr)
library(ggplot2)
library(doMC)
library(rtracklayer)
library(BSgenome)
library('BSgenome.Hsapiens.UCSC.hg18')
library(TxDb.Hsapiens.UCSC.hg18.knownGene)

#### Register rtracklayer backend ####
session <- browserSession()
genome(session) <- "hg18"

#### Define genomic features ####
# The list of chromosomes to be studied - we use chr1, ..., chr22, chrX and chrY. We don't use chrM because it is a circular chromosome, which breaks subsetByOverlaps()
chromosome.list <- list('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                        'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY')
names(chromosome.list) <- chromosome.list

# CpG islands (CGIs)
CGI <- read.table("/home/users/lab0605/hickey/Rafas-cpg-islands-hg18.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
CGI <- GRanges(seqnames = CGI$chr, ranges = IRanges(start = CGI$start, end = CGI$end), 
               CpGcount = CGI$CpGcount, GCcontent = CGI$GCcontent, pctGC = CGI$pctGC, 
               obsExp =  CGI$obsExp)
seqlevels(CGI, force = TRUE) <- names(chromosome.list) # Drop chromosomes not in input data

# Regions outside CGIs
outside.CGI <- gaps(CGI)

# UCSC knownGenes - each range is [transcription start site, transcription end site] and thus includes intronic sequence. Furthermore, there are multiple transcripts included per gene.
query <- ucscTableQuery(session, "knownGene") # generate a query
knownGene <- track(query, asRangedData = FALSE) # which in turn can be used to get the data
seqlevels(knownGene, force = TRUE) <- names(chromosome.list) # Drop chromosomes not in input data

# Regions outside any knownGene
tmp <- knownGene
strand(tmp) <- '*'
seqlengths(tmp) <- NA
outside.knownGene <- gaps(reduce(tmp))
rm(tmp)

# knownGene exons
knownGene.exons <- unlist(exonsBy(TxDb.Hsapiens.UCSC.hg18.knownGene))
seqlevels(knownGene.exons, force = TRUE) <- names(chromosome.list) # Drop chromosomes not in input data

# Regions outside any knownGene exons
tmp <- knownGene.exons
strand(tmp) <- "*"
seqlengths(tmp) <- NA
outside.knownGene.exons <- gaps(reduce(tmp))
rm(tmp)

# Regions inside knownGene but not exonic
nonexonic.knownGene <- subsetByOverlaps(outside.knownGene.exons, knownGene, type = 'within')

# Transcription start sites (TSS)
TSS <- as.vector(ifelse(strand(knownGene) == "+", start(knownGene), end(knownGene)))

# knownGene transcription start sites +/- 500bp
TSS_500 <- GRanges(seqnames = seqnames(knownGene), IRanges(start = pmax(1, TSS - 500), end = TSS + 500))

# knownGene transcription start sites +/- 2000bp
TSS_2000 <- GRanges(seqnames = seqnames(knownGene), IRanges(start = pmax(1, TSS - 2000), end = TSS + 2000))

# Transcription end sites (TES)
TES <- as.vector(ifelse(strand(knownGene) == "+", end(knownGene), start(knownGene)))

# knownGene transcription end sites +/- 500bp
TES_500 <- GRanges(seqnames = seqnames(knownGene), IRanges(start = pmax(1, TES - 500), end = TES + 500))

# knownGene transcription end sites +/- 2000bp
TES_2000 <- GRanges(seqnames = seqnames(knownGene), IRanges(start = pmax(1, TES - 2000), end = TES + 2000))

# CGI shores - defined as 2kb up and downstream of CGIs
CGI.shore <- c(flank(CGI, width = 2000, start = TRUE), flank(CGI, width = 2000, start = FALSE))

#### Remove unneccesary objects and save image ####
rm(chromosome.list, query, session, TES, TSS)
save.image('/home/users/lab0605/hickey/Comethylation_scripts/genomic-features.RData')