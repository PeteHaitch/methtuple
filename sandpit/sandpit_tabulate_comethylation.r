# A script to tabulate the methylation state of reads overlapping two proximal CpGs
# Peter Hickey (hickey@wehi.edu.au)
# 19/12/2011

## This is a sandpit for methods development

### TO DO 
##### Return DataFrame rather than data.frame for all functions
# add an option to chooseCpGPair that favours assigning reads to CpG pairs if they overlap both CpGs in the pair
# Should perhaps be using NA's instead of NULL's for CpG pairs with no comethylation table.
# Change if(){}else{} statements to ifelse() function calls
# Do scatterplots of methylation levels at CpG pairs with lag d
# I think "bookended" CpG pairs are fine - it's the reads that must be independent, not the CpGs that make up the pairs
# Indeed I can probably generalise the whole process to the do the following:
# (0) For every CpG store its "context", e.g. CGI, promoter, gene body, etc.
# (1) For every read select at random two CpGs it overlaps
#   (a) If the read only overlaps one CpG then what to do with it? 
#       Such reads don't contribute to the comethylation odds ratio, but are useful for estimating the marginal methylation probabilities
#   (b) Need to sensibly handle reads that report something other than "Z" or "z" at a CpG
#       May want to chuck these reads out once the tables are created
# (2) Keep a record of all CpG pairs constructed by this process and their corresponding methylation values (store as a hash table?)
# (3) Summarise comethylation by creating tables based on all CpG pairs satisfying some condition(s),
#     e.g. lag = d & both CpGs in a promoter region
# What % of CpGs are covered in the Hansen et al. Illumina data?
## I can do the scatterplot "higher-level" comethylation just based on the pileup.
## We'll later disect to what extent the lower-level comethylation drives the higher-level comethylation.

# 17/12/2011 
##############################################################################################################
library(reshape2)
## Try using the melt() function to reshape the lor list
a <- SRR206931.lor
# I forgot to include the distance variable in the lor lists, so add them now
dist <- as.numeric(rownames(a$chr1)) - 1
b <- melt(a)
b$dist <- rep(dist, 23 * 2)
b <- subset(b, variable == "lor")

# Simple version of what I want
quickBox <- function(lor, sample.name){
  dist <- as.numeric(rownames(lor$chr1)) - 1
  x <- melt(lor)
  x$dist <- rep(dist, 23 * 2)
  x <- subset(x, variable == "lor")
  boxplot(value ~ dist, data = x, xlab = "Distance", ylab = "Log odds ratio", main = sample.name, ylim = c(-4, 7))
}

# Fancy version
slowBox <- function(lor){
  dist <- as.numeric(rownames(lor$chr1)) - 1
  x <- melt(lor)
  x$dist <- rep(dist, 23 * 2)
  x <- subset(x, variable == "lor")
  p <- ggplot(x, aes(factor(dist), value))   
  p + geom_boxplot()
}

# Generate quickBox() plots for all samples
setwd("/usr/local/work/hickey/Hansen_BSseq/sandpit/plots/")
pdf(file = "SRR206931/SRR206931_lor_boxplot.pdf")
quickBox(SRR206931.lor, "SRR206931")
tools::compactPDF("SRR206931/SRR206931_lor_boxplot.pdf")
dev.off()
pdf(file = "SRR206932/SRR206932_lor_boxplot.pdf")
quickBox(SRR206932.lor, "SRR206932")
tools::compactPDF("SRR206932/SRR206932_lor_boxplot.pdf")
dev.off()
pdf(file = "SRR206933/SRR206933_lor_boxplot.pdf")
quickBox(SRR206933.lor, "SRR206933")
tools::compactPDF("SRR206933/SRR206933_lor_boxplot.pdf")
dev.off()
pdf(file = "SRR206934/SRR206934_lor_boxplot.pdf")
quickBox(SRR206934.lor, "SRR206934")
tools::compactPDF("SRR206934/SRR206934_lor_boxplot.pdf")
dev.off()
pdf(file = "SRR206935/SRR206935_lor_boxplot.pdf")
quickBox(SRR206935.lor, "SRR206935")
tools::compactPDF("SRR206935/SRR206935_lor_boxplot.pdf")
dev.off()
pdf(file = "SRR207940/SRR207940_lor_boxplot.pdf")
quickBox(SRR207940.lor, "SRR207940")
tools::compactPDF("SRR207940/SRR207940lor_boxplot.pdf")
dev.off()

# 14/12/2011 PM
##############################################################################################################
# Boxplots of log odds ratios - one boxplot per "distance/lag" and each chromosome contributing one datapoint to each boxplot,
# i.e. each boxplot is built on 22 datapoints (one for each of chr1,...,chr22 since I am not currently applying it to chrX or other chromosomes
source("tabulate_comethylation.R")
load("SRR206931_comethylation.RData")
library(ggplot2)

# Take methlyation data.frame and turn into a comethylation table
methylDF2Table <- function(x, type = "collapsed"){
  cometh.data.list <- values(x)$cometh.data
  cometh.data.joined <- Reduce(function(y, z) rbind(y, z), cometh.data.list, accumulate=F)
  if(type == "all"){
    # Create a table for every CpG pair in cometh.data.joined i.e. without collapsing
    table.all <- table(cometh.data.joined)
    return(table.all)
  }
  # Create a table by collapsing over all CpG pairs in cometh.data.joined
  table.collapsed <- table(cometh.data.joined[, 1:2])
  if(type == "collapsed"){
    return(table.collapsed)
  } else if(type == "subtable"){
    # Create a sub-table of table.tall. The sub-table is conveniently ordered via the [-function
    subtable <- table.collapsed[c("Z", "z"), c("Z", "z")]
  return(subtable)
  }
}

# Checks to make sure entry exists, if it doesn't assigns a value of 0.5 = 0 + 0.5
logOdds <- function(x){
  aa <- ifelse("Z" %in% rownames(x) & "Z" %in% colnames(x), x["Z", "Z"] + 0.5, 0 + 0.5)
  bb <- ifelse("Z" %in% rownames(x) & "z" %in% colnames(x), x["Z", "z"] + 0.5, 0 + 0.5)
  cc <- ifelse("z" %in% rownames(x) & "Z" %in% colnames(x), x["z", "Z"] + 0.5, 0 + 0.5)
  dd <- ifelse("z" %in% rownames(x) & "z" %in% colnames(x), x["z", "z"] + 0.5, 0 + 0.5)
  lor <- log((aa*dd)/(bb*cc))
  ASE <- sqrt(1/aa + 1/bb + 1/cc + 1/dd) # See pp54 of Agresti
  return(c(lor, ASE))
}

## TODO: Separate the computation of the tables from the plotting into two distinct functions
# A function to compute the log odds ratio for a single chromosome and (optionally) plot the result.
calcORs <- function(chr.name, sample, sample.name, plot = FALSE){
  chr.data <- sample[[chr.name]]
  chr.data <- split(chr.data, width(chr.data))
  lst <- as(chr.data, "list")
  # Probably don't even need 'subtable' option since I match on column- and row-names, anyway.
  tables <- mclapply(X = lst, FUN = methylDF2Table, type = "collapsed")
  # Compute the log-odds-ratio for a range of distances
  n <- 80 # Max distance 
  df <- data.frame(lor = rep(NA, n), ASE = rep(NA, n))
  for(i in 1:n){
    df[i,] <- logOdds(tables[[i]])
  }
  df$distance <- as.numeric(names(tables)[1:80]) - 4
  if(plot == TRUE){
    limits <- aes(ymax = lor + 2*ASE, ymin = lor - 2*ASE)
    filename <- paste(sample.name, "_", chr.name, ".pdf", sep = "")
    m <- ggplot(data = df, aes(x = distance, y = lor)) + geom_point() + geom_errorbar(limits) + opts(title = paste(sample.name, ": ", chr.name, sep = ""))
    ggsave(filename = filename, plot = m)
    dev.off()
  }
  return(df)
}

# Apply to all chromosomes for SRR206931 but don't produce any plots
setwd("/usr/local/work/hickey/Hansen_BSseq/sandpit/plots/SRR206931")
chr.names <- vector("list", 23)
for(i in 1:22){
  chr.names[[i]] <- paste("chr", i, sep = "")
}
chr.names[[23]] <- "chrX"

SRR206931.lor <- lapply(X=chr.names, FUN = calcORs, sample = SRR206931.tab, sample.name = "SRR206931", plot = FALSE) # Currently running
names(SRR206931.lor) <- chr.names
save(SRR206931.lor, file = "../../SRR206931_lor.RData")

# Repeat for SRR206932, etc.
load("/usr/local/work/hickey/Hansen_BSseq/sandpit/SRR206932_comethylation.RData")
setwd("/usr/local/work/hickey/Hansen_BSseq/sandpit/plots/SRR206932")
SRR206932.lor <- lapply(X=chr.names, FUN = calcORs, sample = SRR206932.tab, sample.name = "SRR206932", plot = FALSE) # Currently running
names(SRR206932.lor) <- chr.names
save(SRR206932.lor, file = "../../SRR206932_lor.RData")

load("/usr/local/work/hickey/Hansen_BSseq/sandpit/SRR206933_comethylation.RData")
setwd("/usr/local/work/hickey/Hansen_BSseq/sandpit/plots/SRR206933")
SRR206933.lor <- lapply(X=chr.names, FUN = calcORs, sample = SRR206933.tab, sample.name = "SRR206933", plot = FALSE) # Currently running
names(SRR206933.lor) <- chr.names
save(SRR206933.lor, file = "../../SRR206933_lor.RData")

load("/usr/local/work/hickey/Hansen_BSseq/sandpit/SRR206934_comethylation.RData")
setwd("/usr/local/work/hickey/Hansen_BSseq/sandpit/plots/SRR206934")
SRR206934.lor <- lapply(X=chr.names, FUN = calcORs, sample = SRR206934.tab, sample.name = "SRR206934", plot = FALSE) # Currently running
names(SRR206934.lor) <- chr.names 
save(SRR206934.lor, file = "../../SRR206934_lor.RData")

load("/usr/local/work/hickey/Hansen_BSseq/sandpit/SRR206935_comethylation.RData")
setwd("/usr/local/work/hickey/Hansen_BSseq/sandpit/plots/SRR206935")
SRR206935.lor <- lapply(X=chr.names, FUN = calcORs, sample = SRR206935.tab, sample.name = "SRR206935", plot = FALSE) # Currently running
names(SRR206935.lor) <- chr.names 
save(SRR206935.lor, file = "../../SRR206935_lor.RData")

load("/usr/local/work/hickey/Hansen_BSseq/sandpit/SRR207940_comethylation.RData")
setwd("/usr/local/work/hickey/Hansen_BSseq/sandpit/plots/SRR207940")
SRR207940.lor <- lapply(X=chr.names, FUN = calcORs, sample = SRR207940.tab, sample.name = "SRR207940", plot = FALSE) # Currently running
names(SRR207940.lor) <- chr.names 
save(SRR207940.lor, file = "../../SRR207940_lor.RData") # UP TO HERE



# 14/12/2011 AM
##############################################################################################################
source("tabulate_comethylation.R")
load("SRR206932_comethylation.RData")
library(ggplot2)

# Take methlyation data.frame and turn into a comethylation table
methylDF2Table <- function(x, type = "collapsed"){
  cometh.data.list <- values(x)$cometh.data
  cometh.data.joined <- Reduce(function(y, z) rbind(y, z), cometh.data.list, accumulate=F)
  if(type == "all"){
    # Create a table for every CpG pair in cometh.data.joined i.e. without collapsing
    table.all <- table(cometh.data.joined)
    return(table.all)
  }
  # Create a table by collapsing over all CpG pairs in cometh.data.joined
  table.collapsed <- table(cometh.data.joined[, 1:2])
  if(type == "collapsed"){
    return(table.collapsed)
  } else if(type == "subtable"){
    # Create a sub-table of table.tall. The sub-table is conveniently ordered via the [-function
    subtable <- table.collapsed[c("Z", "z"), c("Z", "z")]
  return(subtable)
  }
}

logOdds <- function(x){
  aa <- x["Z", "Z"] + 0.5
  bb <- x["Z", "z"] + 0.5
  cc <- x["z", "Z"] + 0.5
  dd <- x["z", "z"] + 0.5
  lor <- log((aa*dd)/(bb*cc))
  ASE <- sqrt(1/aa + 1/bb + 1/cc + 1/dd) # See pp54 of Agresti
  return(c(lor, ASE))
}

# A function to plot a single chromosome's data
plotORs <- function(chr.name, sample, sample.name){
  chr.data <- sample[[chr.name]]
  chr.data <- split(chr.data, width(chr.data))
  lst <- as(chr.data, "list")
  tables <- mclapply(X = lst, FUN = methylDF2Table, type = "collapsed")
  # Compute the log-odds-ratio for a range of distances
  n <- 80 # Max distance 
  df <- data.frame(lor = rep(NA, n), ASE = rep(NA, n))
  for(i in 1:n){
    df[i,] <- logOdds(tables[[i]])
  }
  df$distance <- as.numeric(names(tables)[1:80]) - 4
  limits <- aes(ymax = lor + 2*ASE, ymin = lor - 2*ASE)
  filename <- paste(sample.name, "_", chr.name, ".pdf", sep = "")
  m <- ggplot(data = df, aes(x = distance, y = lor)) + geom_point() + geom_errorbar(limits) + opts(title = paste(sample.name, ": ", chr.name, sep = ""))
  ggsave(filename = filename, plot = m)
  dev.off()
}

# Apply to all chromosomes
setwd("/usr/local/work/hickey/Hansen_BSseq/sandpit/plots/SRR206932")
chr.names <- vector("list", 22)
for(i in 1:22){
  chr.names[[i]] <- paste("chr", i, sep = "")
}

lapply(X=chr.names, FUN = plotORs, sample = SRR206932.tab, sample.name = "SRR206932")


# 13/12/2011
##############################################################################################################

##### NEED TO MAKE THE AGGREGATION STEP MORE ROBUST, e.g. it will break if there are tables that don't contain Z's and z's

## Produce graphs of odds ratio +/- two standard errors for each chromsome and sample, stratified by pair-distance
## Basically try to turn the work on 06/12/2011 (AM) into a function
source("tabulate_comethylation.R")
load("SRR206931_comethylation.RData")
library(ggplot2)

# Take methlyation data.frame and turn into a comethylation table
methylDF2Table <- function(x, type = "collapsed"){
  cometh.data.list <- values(x)$cometh.data
  cometh.data.joined <- Reduce(function(y, z) rbind(y, z), cometh.data.list, accumulate=F)
  if(type == "all"){
    # Create a table for every CpG pair in cometh.data.joined i.e. without collapsing
    table.all <- table(cometh.data.joined)
    return(table.all)
  }
  # Create a table by collapsing over all CpG pairs in cometh.data.joined
  table.collapsed <- table(cometh.data.joined[, 1:2])
  if(type == "collapsed"){
    return(table.collapsed)
  } else if(type == "subtable"){
    # Create a sub-table of table.tall. The sub-table is conveniently ordered via the [-function
    subtable <- table.collapsed[c("Z", "z"), c("Z", "z")]
  return(subtable)
  }
}

## A function to compute the log-odds-ratio and its asymptotic standard error
## Basically use the sub-table that only involved (Z, z)'s = (M, U)'s and compute the log odds ratio 
## Compute the modified unconditional maximum likelihood estimator of lambda (the log odds ratio)
## See pp68 of Plackett
logOdds <- function(x){
  aa <- x["Z", "Z"] + 0.5
  bb <- x["Z", "z"] + 0.5
  cc <- x["z", "Z"] + 0.5
  dd <- x["z", "z"] + 0.5
  lor <- log((aa*dd)/(bb*cc))
  ASE <- sqrt(1/aa + 1/bb + 1/cc + 1/dd) # See pp54 of Agresti
  return(c(lor, ASE))
}

# A function to plot a single chromosome's data
plotORs <- function(chr.name, sample, sample.name){
  chr.data <- sample[[chr.name]]
  chr.data <- split(chr.data, width(chr.data))
  lst <- as(chr.data, "list")
  tables <- mclapply(X = lst, FUN = methylDF2Table, type = "collapsed")
  # Compute the log-odds-ratio for a range of distances
  n <- 80 # Max distance 
  df <- data.frame(lor = rep(NA, n), ASE = rep(NA, n))
  for(i in 1:n){
    df[i,] <- logOdds(tables[[i]])
  }
  df$distance <- as.numeric(names(tables)[1:80]) - 4
  limits <- aes(ymax = lor + 2*ASE, ymin = lor - 2*ASE)
  filename <- paste(sample.name, "_", chr.name, ".pdf", sep = "")
  m <- ggplot(data = df, aes(x = distance, y = lor)) + geom_point() + geom_errorbar(limits) + opts(title = paste(sample.name, ": ", chr.name, sep = ""))
  ggsave(filename = filename, plot = m)
  dev.off()
}

# A test of the plotORs function
plotORs("chr22", SRR206931.tab, "SRR206931") # Seems to work

# Apply to all chromosomes
setwd("/usr/local/work/hickey/Hansen_BSseq/sandpit/plots/")
chr.names <- vector("list", 23)
for(i in 1:22){
  chr.names[[i]] <- paste("chr", i, sep = "")
}

# chr21 doesn't work - figure out why
lapply(X=chr.names, FUN = plotORs, sample = SRR206931.tab, sample.name = "SRR206931")


# 06/12/2011 PM
##############################################################################################################
# Let's process some more samples
source("tabulate_comethylation.R")
options(cores = 20)

hg19.cpgs <- import("/usr/local/work/hickey/CpGs/hg19_CpGs.bed.gz")
cpg.pairs <- makeCpGPairs(cpgs = hg19.cpgs, type = "disjoint")

SRR206932.reads <- readBismarkBam("/usr/local/work/hickey/Hansen_BSseq/bams/Illumina_bismark/SRR206932_bismark_sorted.bam")
SRR206932.tab <- makeCoMethTable(reads = SRR206932.reads, cpg.pairs =  cpg.pairs, pair.choice = "random")
rm(SRR206932.reads)
save(SRR206931.tab, file = "SRR206932_comethylation.RData")

SRR206933.reads <- readBismarkBam("/usr/local/work/hickey/Hansen_BSseq/bams/Illumina_bismark/SRR206933_bismark_sorted.bam")
SRR206933.tab <- makeCoMethTable(reads = SRR206933.reads, cpg.pairs =  cpg.pairs, pair.choice = "random")
rm(SRR206933.reads)
save(SRR206933.tab, file = "SRR206933_comethylation.RData")

SRR206934.reads <- readBismarkBam("/usr/local/work/hickey/Hansen_BSseq/bams/Illumina_bismark/SRR206934_bismark_sorted.bam")
SRR206934.tab <- makeCoMethTable(reads = SRR206934.reads, cpg.pairs =  cpg.pairs, pair.choice = "random")
rm(SRR206934.reads)
save(SRR206934.tab, file = "SRR206934_comethylation.RData")

SRR206935.reads <- readBismarkBam("/usr/local/work/hickey/Hansen_BSseq/bams/Illumina_bismark/SRR206935_bismark_sorted.bam")
SRR206935.tab <- makeCoMethTable(reads = SRR206935.reads, cpg.pairs =  cpg.pairs, pair.choice = "random")
rm(SRR206935.reads)
save(SRR206935.tab, file = "SRR206935_comethylation.RData")

SRR207940.reads <- readBismarkBam("/usr/local/work/hickey/Hansen_BSseq/bams/Illumina_bismark/SRR207940_bismark_sorted.bam")
SRR207940.tab <- makeCoMethTable(reads = SRR207940.reads, cpg.pairs =  cpg.pairs, pair.choice = "random")
rm(SRR207940.reads)
save(SRR207940.tab, file = "SRR207940_comethylation.RData")

# 06/12/2011 AM
##############################################################################################################
source("tabulate_comethylation.R")
load("SRR206931_comethylation.RData")
## Merge the methylation calls by CpG pair "internal-distance"
# Grab a chromosome to play with
a <- SRR206931.tab$chr13
## Create tables for every CpG pair
#a.tab <- mclapply(elementMetadata(a)$cometh.data,  table)
#elementMetadata(a)$cometh.tab <- a.tab
# Split the data by "internal-distance"
b <- split(a, width(a))
# How many splits are there?
length(names(b))

# rbind all methylation calls together per "internal-distance"
d <- b$"4"

#### PROBLEM - missing in both reads! 
#### REASON - Reads mapped to the Crick strand can overlap the G but not the C in CpG_1 and then not overlap at all CpG_2 in the pair.
#### Schematic (Crick strand)
#### 5'.....CG....CG......3' (reference strand is Watson strand)
####  <-----|                (an example of the read described above: reports "-", "-" for this pair of CpGs)
#### Conversely, reads mapped to the Watson strand overlap the G but not the C in CpG_2 and then not overlap at all CpG_1 in the pair.
#### Schematic (Watson strand)
#### 5'.....CG....CG......3' (reference strand is Watson strand)
####               |------> (an example of the read described above: reports "-", "-" for this pair of CpGs)
values(d)$cometh.data[[2237]] # chr13: 34185374-34185377

e <- values(d)$cometh.data

tmp <- Reduce(function(x, y) rbind(x, y), e, accumulate=F)
# Create a table by collapsing over all CpG pairs with width = 4
tmp.table <- table(tmp[, 1:2])
# Create a table for every CpG pair with width = 4 i.e. without collapsing
tmp.table.all <- table(tmp)
# Create a sub-table. The sub-table can be conveniently ordered via the [-function
tmp2 <- tmp.table[c("Z", "z"), c("Z", "z")]

## Turn the above into a function
methylDF2Table <- function(x, type = "collapsed"){
  cometh.data.list <- values(x)$cometh.data
  cometh.data.joined <- Reduce(function(y, z) rbind(y, z), cometh.data.list, accumulate=F)
  if(type == "all"){
    # Create a table for every CpG pair in cometh.data.joined i.e. without collapsing
    table.all <- table(cometh.data.joined)
    return(table.all)
  }
  # Create a table by collapsing over all CpG pairs in cometh.data.joined
  table.collapsed <- table(cometh.data.joined[, 1:2])
  if(type == "collapsed"){
    return(table.collapsed)
  } else if(type == "subtable"){
    # Create a sub-table of table.tall. The sub-table is conveniently ordered via the [-function
    subtable <- table.collapsed[c("Z", "z"), c("Z", "z")]
  return(subtable)
  }
}
#methylDF2Table(b$"4", type = "all")
#methylDF2Table(b$"4", type = "collapsed")
#methylDF2Table(b$"4", type = "subtable")

## Let's do some time-comparisons using mclapply to see if it's worth the effort of converting the GRangesList to a list
## Looks like it is worth doing even when only using 10 cores
tmp <- lapply(X = b, FUN = methylDF2Table, type = "collapsed")
# See how long it takes to process chromosome 13 for SRR206931
system.time(lapply(X = b, FUN = methylDF2Table, type = "collapsed"))
#user  system elapsed 
#734.833   0.621 736.111
system.time(as(b, "list"))
#user  system elapsed 
#165.796   1.542 167.716 
lst <- as(b, "list")
options(cores = 24)
tmp2 <- mclapply(X = lst, FUN = methylDF2Table, type = "collapsed")
system.time(mclapply(X = lst, FUN = methylDF2Table, type = "collapsed"))
#user  system elapsed 
#620.996  18.471  44.872
options(cores = 10)
system.time(mclapply(X = lst, FUN = methylDF2Table, type = "collapsed"))
#user  system elapsed 
#513.000   7.259  80.359 
options(cores = 40)
system.time(mclapply(X = lst, FUN = methylDF2Table, type = "collapsed"))
# user  system elapsed 
#684.320  31.702  35.425 
options(cores = 10)
tables <- mclapply(X = lst, FUN = methylDF2Table, type = "collapsed")

# It looks like only (Z, z, ., "") are possible values but not all tables are 4x4
row.names <- table(unlist(lapply(tables[1:length(tables)], rownames))) # What row-types are there?
col.names <- table(unlist(lapply(tables[1:length(tables)], function(x){rownames(t(x))}))) # What column-types are there?
#chisq.test(rbind(row.names, col.names))
dims <- lapply(tables, dim) # Not all tables are 4x4

# A table for testing inferential procedures
test.tab <- tables[[1]]
row.names(test.tab) <- c("Z", "z", ".", "-")
colnames(test.tab) <- c("Z", "z", "-", ".")
test.tab <- test.tab[c("Z", "z", ".", "-"), c("Z", "z", ".", "-")]

sub.test.tab <- test.tab[c("Z", "z"), c("Z", "z")]
df <- data.frame(CpG1 = c("Z", "Z", "z", "z"), CpG2 = c("Z", "z", "Z", "z"), 
                 counts = c(test.tab["Z", "Z"], test.tab["Z", "z"], test.tab["z", "Z"], test.tab["z", "z"]))
glm.1 <- glm(counts ~ CpG1 + CpG2, family = poisson, data = df)
glm.2 <- glm(counts ~ CpG1 * CpG2, family = poisson, data = df)
glm.3 <- glm(Freq ~ first + second, family = poisson, data = as.data.frame(test.tab))
glm.4 <- glm(Freq ~ first * second, family = poisson, data = as.data.frame(test.tab))

## See what happens when I add in the all reads that don't overlap either CpG_1 or CpG_2 
## There 2888270 reads for sample SRR20691, minus the sum(tt) that overlap this pair of CpGs
tt <- test.tab
tt[4, 3] <- 2888270 - sum(tt) + tt[4,3]
glm.5 <- glm(Freq ~ first + second, family = poisson, data = as.data.frame(tt))
glm.6 <- glm(Freq ~ first * second, family = poisson, data = as.data.frame(tt))

## Basically use the sub-table that only involved (Z, z)'s = (M, U)'s and compute the log odds ratio 
## Compute the modified unconditional maximum likelihood estimator of lambda (the log odds ratio)
## See pp68 of Plackett
logOdds <- function(x){
  aa <- x["Z", "Z"] + 0.5
  bb <- x["Z", "z"] + 0.5
  cc <- x["z", "Z"] + 0.5
  dd <- x["z", "z"] + 0.5
  lor <- log((aa*dd)/(bb*cc))
  ASE <- sqrt(1/aa + 1/bb + 1/cc + 1/dd) # See pp54 of Agresti
  return(c(lor, ASE))
}

# Compute this for a bunch of distances
n <- 80
x <- data.frame(lor = rep(NA, n), ASE = rep(NA, n))
for(i in 1:n){
  x[i,] <- logOdds(tables[[i]])
}

# Plot this
library(ggplot2)
df <- x
df$distance <- names(tables)[1:80]
limits <- aes(ymax = lor + 2*ASE, ymin = lor - 2*ASE)
qplot(x = distance, y = lor, data = df) + geom_errorbar(limits)  # Add a smoothed line
ggplot(data = df, aes(x = distance, y = lor)) + geom_point() + geom_errorbar(limits) + stat_smooth(method=loess, se = FALSE, span = 0.5)

# Or to plot it on the back-transformed scale
limits <- aes(ymax = exp(lor + 2*ASE), ymin = exp(lor - 2*ASE))
qplot(y = exp(lor), data = df) + geom_errorbar(limits) # Add a smoothed line

## Do I want use the conditional or unconditional odds ratio estimator?
## Currently using the continuity corrected sample odds ratio
# Compute the continuity corrected odds ratio and it's standard error for a 2x2 table
OR <- function(x, a, b, r, s){
  aa <- x[a,b] + 0.5
  bb <- x[a,s] + 0.5
  cc <- x[r,b] + 0.5
  dd <- x[r,s] + 0.5
  or <- (aa*dd)/(bb*cc)
  se <- 1/aa + 1/bb + 1/cc + 1/dd
  #return(data.frame(or = or, se = se))
  return(c(or, se))
}

# log odds ratio
LOR <- function(x, a, b, r, s){
  aa <- x[a,b] + 0.5
  bb <- x[a,s] + 0.5
  cc <- x[r,b] + 0.5
  dd <- x[r,s] + 0.5
  lor <- log((aa*dd)/(bb*cc))
  lse <- log(1/aa + 1/bb + 1/cc + 1/dd)
  return(c(lor, lse))
}

lambda <- c(LOR(test.tab, 1, 1, 4, 4)[1], LOR(test.tab, 2, 1, 4, 4)[1], LOR(test.tab, 3, 1, 4, 4)[1], 
            LOR(test.tab, 1, 2, 4, 4)[1], LOR(test.tab, 2, 2, 4, 4)[1], LOR(test.tab, 3, 2, 4, 4)[1], 
            LOR(test.tab, 1, 3, 4, 4)[1], LOR(test.tab, 2, 3, 4, 4)[1], LOR(test.tab, 3, 3, 4, 4)[1])
lambda.matrix <- matrix(lambda, nrow = 3, byrow = F)


# 05/12/2011
##############################################################################################################
## Merge the methylation calls by CpG pair "internal-distance"
# Grab a chromosome to play with
a <- SRR206931.tab$chr13
## Create tables for every CpG pair
#a.tab <- mclapply(elementMetadata(a)$cometh.data,  table)
#elementMetadata(a)$cometh.tab <- a.tab
# Split the data by "internal-distance"
b <- split(a, width(a))
# How many splits are there?
length(names(b))

# rbind all methylation calls together per "internal-distance"
d <- b$"4"
e <- values(d)$cometh.data




# 30/11/2011 PM
# binfbig1 code
### THIS WORKS
##############################################################################################################
source("tabulate_comethylation.R")


hg19.cpgs <- import("/usr/local/work/hickey/CpGs/hg19_CpGs.bed.gz")
cpg.pairs <- makeCpGPairs(cpgs = hg19.cpgs, type = "disjoint")

SRR206931.reads <- readBismarkBam("/usr/local/work/hickey/Hansen_BSseq/bams/Illumina_bismark/SRR206931_bismark_sorted.bam")
SRR206931.tab <- makeCoMethTable(reads = SRR206931.reads, cpg.pairs =  cpg.pairs, pair.choice = "random")
rm(SRR206931.reads)

# 29/11/2011 PM
# Local code
##############################################################################################################
setwd("~/Desktop/Comethylation/sandpit/")
source("~/Desktop/Comethylation/R/tabulate_comethylation.R")

### Test data set
# Read in the bam file
reads <- readBismarkBam("SRR206931_bismark_sorted.bam")

# Check if there is a genome-wide strand bias
table(reads@strand)/length(reads)  * 100

# Import all chr13 CpGs from a bed file
bed <- import("hg19_CpGs_chr13.bed")

# Convert CpGs to CpG pairs
cpg.pairs <- makeCpGPairs(bed)
cpg.pairs <- cpg.pairs$chr13

# Test the chooseCpGPair function
a <- chooseCpGPair(reads = reads, cpg.pairs = cpg.pairs, type = "random")

# Test the olaps2List function
b <- olaps2List(a)
d <- b[1:10]

chr.names <- "chr13"

# makeTable2 can only handle a single (named) chromosome
makeTable2 <- function(chr, reads, cpg.pairs, input.for.makeTable){ 
  chr.data <- input.for.makeTable
  cpg.pairs <- cpg.pairs
  mclapply(X = names(chr.data), FUN = gmc, overlaps.list = chr.data, reads = reads, cpg.pairs = cpg.pairs)
}

# A list-ified version of getMethylCallOld()
gmc <- function(name, overlaps.list, reads, cpg.pairs){
  reads <- reads[as.vector(unlist(overlaps.list[name]))]
  cpg.pair <- cpg.pairs[as.numeric(name), ]
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
  # Combine the results and return as a data.frame 
  # Need to name each data.frame by a CpG pair identifier
  # This is a kludge that stores index of the cpg.pair in the first column name
   output <- rbind(watson, crick)
   names(output) <- name
   return(output)
  #return(rbind(watson, crick))
}

e <- makeTable2(chr=chr.names, reads=reads, cpg.pairs=cpg.pairs, input.for.makeTable=d)
# f <- makeTable2(chr=chr.names, reads=reads, cpg.pairs=cpg.pairs, input.for.makeTable=b)

# Now add comethylation results to the cpg.pairs object
cpg.pairs.gr <- as(cpg.pairs, "GRanges")
# A function to extract the index for the cpg.pair from the comethylation table
getID <- function(x){
  as.numeric(names(x)[1])
}

IDs <- unlist(lapply(f, getID))
cpg.pairs.gr@elementMetadata$cometh.data <- NA
cpg.pairs.gr@elementMetadata$cometh.data[IDs] <- f 




# 29/11/2011 AM
# Code for binfbig1
##############################################################################################################
source("tabulate_comethylation.R")

hg19.cpgs <- import("/usr/local/work/hickey/CpGs/hg19_CpGs.bed.gz")
cpg.pairs <- makeCpGPairs(cpgs = hg19.cpgs, type = "disjoint")

SRR206931.reads <- readBismarkBam("/usr/local/work/hickey/Hansen_BSseq/bams/Illumina_bismark/SRR206931_bismark_sorted.bam")
SRR206931.tab <- makeCoMethTable(reads = SRR206931.reads, cpg.pairs =  cpg.pairs, pair.choice = "random")
rm(SRR206931.reads)

SRR206932.reads <- readBismarkBam("/usr/local/work/hickey/Hansen_BSseq/bams/Illumina_bismark/SRR206932_bismark_sorted.bam")
SRR206932.tab <- makeCoMethTable(reads = SRR206932.reads, cpg.pairs =  cpg.pairs, pair.choice = "random")
rm(SRR206932.reads)

SRR206933.reads <- readBismarkBam("/usr/local/work/hickey/Hansen_BSseq/bams/Illumina_bismark/SRR206933_bismark_sorted.bam")
SRR206933.tab <- makeCoMethTable(reads = SRR206933.reads, cpg.pairs =  cpg.pairs, pair.choice = "random")
rm(SRR206933.reads)

SRR206934.reads <- readBismarkBam("/usr/local/work/hickey/Hansen_BSseq/bams/Illumina_bismark/SRR206934_bismark_sorted.bam")
SRR206934.tab <- makeCoMethTable(reads = SRR206934.reads, cpg.pairs =  cpg.pairs, pair.choice = "random")
rm(SRR206934.reads)

SRR206935.reads <- readBismarkBam("/usr/local/work/hickey/Hansen_BSseq/bams/Illumina_bismark/SRR206935_bismark_sorted.bam")
SRR206935.tab <- makeCoMethTable(reads = SRR206935.reads, cpg.pairs =  cpg.pairs, pair.choice = "random")
rm(SRR206935.reads)

SRR207940.reads <- readBismarkBam("/usr/local/work/hickey/Hansen_BSseq/bams/Illumina_bismark/SRR207940_bismark_sorted.bam")
SRR207940.tab <- makeCoMethTable(reads = SRR207940.reads, cpg.pairs =  cpg.pairs, pair.choice = "random")
rm(SRR207940.reads)
# 26/11/2011- 28/11/2011
##############################################################################################################
### To be run on binfbig1
# Let's try to do this genome-wide
### To be run on binfbig1
source("tabulate_comethylation.R")
reads <- readBismarkBam("/usr/local/work/hickey/Hansen_BSseq/bams/Illumina_bismark/SRR206931_bismark_sorted.bam")
bed <- import("/usr/local/work/hickey/CpGs/hg19_CpGs.bed.gz")
bed.split <- split(bed, space(bed), drop = TRUE) # Only used for testing

# A chromosome with an even number of CpGs
even <- bed.split[[6]]
even.pairs <- makeCpGPairs(even)
even.pairs2 <- makeCpGPairs2(even)
# A chromosome with an odd number of CpGs
odd <- bed.split[[3]]
odd.pairs <- makeCpGPairs(odd)
odd.pairs2 <- makeCpGPairs2(odd)
## An even number of CpGs makes an odd number of CpG pairs (using makeCpGPairs())
## An odd number of CpGs makes an even number of CpG pairs (using makeCpGPairs())

# Make bed into a cpg.pairs object
params.cpg.pairs <- RDApplyParams(rangedData = bed, applyFun = makeCpGPairs)
genome.cpg.pairs <- rdapply(x = params.cpg.pairs)

# Find overlaps of genome.cpg.pairs and reads and choose a single pair of CpG overlaps
random.pairs.list <- lapply(genome.cpg.pairs, chooseCpGPair, reads = reads, type = "random")
leftmost.pairs.list <- lapply(genome.cpg.pairs, chooseCpGPair, reads = reads, type = "first")

# Apply olaps2List() to each element of the output of chooseCpGPair()
input.for.makeTable <- lapply(random.pairs.list, olaps2List) # Can probably mclapply this but may not be worth it

tmp <- input.for.makeTable$chr13

# A function to process a single chromosome
### NEED TO ADD input.for.makeTable as an argument for makeTable
makeTable <- function(chr, reads, all.cpg.pairs){ 
  chr.data <- input.for.makeTable[[chr]]
  cpg.pairs <- all.cpg.pairs[[chr]]
  mclapply(X = names(chr.data), FUN = gmc2, olapsList = chr.data, reads = reads, cpg.pairs = cpg.pairs)
}

z <- makeTable(chr = "chr13", reads = reads, all.cpg.pairs = genome.cpg.pairs)

# Then lapply() makeTable() to two chromosomes
tmp <- input.for.makeTable[1:2]
zz <- lapply(names(tmp), makeTable, reads = reads, all.cpg.pairs = genome.cpg.pairs)

## Apply to whole genome!
### Write all the above as a function

SRR206932 <- lapply(names(input.for.makeTable), makeTable, reads = reads, all.cpg.pairs = genome.cpg.pairs)

# 25/11/2011 PM
##############################################################################################################
# Want to modify getMethylCall() to take in pre-subsetted reads from an olaps2List() instance
# Probably need to pass cpg.pairs and reads to this function as well
gmc2 <- function(name, olapsList, reads, cpg.pairs){
  reads <- reads[as.vector(unlist(olapsList[name]))]
  cpg.pair <- cpg.pairs[as.numeric(name), ]
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
  # Combine the results and return as a data.frame 
  return(rbind(watson, crick))
}

f <- d[1:500]
lapply(X = names(f), FUN = gmc2, f)

mclapply(X = names(f), FUN = gmc2, f)
test <- mclapply(X = names(d), FUN = gmc2, d)

##############################################################################################################

# 25/11/2011 AM
##############################################################################################################
### Test data set
# Read in the bam file
reads <- readBismarkBam("SRR206931_bismark_sorted.bam")

# Check if there is a genome-wide strand bias
table(reads@strand)/length(reads)  * 100

# Import all chr13 CpGs from a bed file
bed <- import("hg19_CpGs_chr13.bed")

# Convert CpGs to CpG pairs
cpg.pairs <- makeCpGPairs(bed)

# Test the chooseCpGPair function
a <- chooseCpGPair(reads = reads, cpg.pairs = cpg.pairs, type = "random")
b <- chooseCpGPair(reads = reads, cpg.pairs = cpg.pairs, type = "first")

# Test the olaps2List function
d <- olaps2List(a)
e <- olaps2List(b)

# How to find out which CpG pair the reads in an element of olaps2List() object overlap - CLUNKY
cpg.pairs[as.numeric(names(d)[5]), ]

## Probably want to randomly assign reads to CpG pairs if they overlap more than one pair, else it will be biased
## e.g.chr13:34,115,480-34,115,626 on cpg.pairs[163536,]
## No reads for CpG_1 but > 1000 for CpG_2 means those 1000 reads couldn't be used for the next CpG if the read overlapped the next CpG too.

# Check out multicore::mcapply()

##############################################################################################################


##############################################################################################################
# 23/11/2011 Problem fixed by Felix and BAM files regenerated using new version of bismark2SAM.pl
# 22/11/2011 Found problem in converted from Bismark to BAM
# Reads are being reverse-complemented in the conversion but the accompanying Phred quality scores and XM tags
# are not being similarly reversed. This breaks my pipeline as I rely on the strand information in the bam file
# to process the XM tag.
# Emailed Felix (the author of Bismark) on 22/11/2011 to ask him to patch the conversion script. Awaiting reply
##############################################################################################################

