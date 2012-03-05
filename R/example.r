# Example usage of tabulate_comethylation.r
# Source the script
source("tabulate_comethylation.R")
# Set the number of cores
options(cores = 20)

# Load the reference genome CpGs
hg19.cpgs <- import("/usr/local/work/hickey/CpGs/hg19_CpGs.bed.gz")
# Construct the CpG pairs object
cpg.pairs <- makeCpGPairs(cpgs = hg19.cpgs, type = "overlapping")
# Calculate the width of the CpG pairs (i.e. the distances between CpG_1 and CpG_2 in each pair)
cpg.pairs.widths <- unlist(lapply(X = cpg.pairs, FUN = width))

# See how many CpG pairs have an intra-pair distance less the upper bound for investigating comethylation
max.read.length <- 87 # Length of the longest read in the library

# keepPair() returns TRUE if the CpG-pair is to be retained for the comethlyation analysis
keepPair <- function(x, read.length){
  keepers <- width(x) <= (read.length - 2)
  x[keepers, ]
}
# Subset the CpG-pairs object by excluding those pairs with intra-pair distance > the upper bound
cpg.pairs.close <- mclapply(X = cpg.pairs, FUN = keepPair, read.length = max.read.length)
# What % of pairs are within the upper-bound
Reduce("+", lapply(cpg.pairs.close, nrow)) / Reduce("+", lapply(cpg.pairs, nrow)) * 100

# Read in the bismark BAM file
example.reads <- readBismarkBam("~/Lister_2009_BS-seq_data/IMR90/post-trim/FASTQ/SRX006783/example_sorted.bam")
example.tab <- makeCoMethTable(reads = example.reads, cpg.pairs =  cpg.pairs, pair.choice = "random")
rm(example.reads)
save(example.tab, file = "example_comethylation.RData")

