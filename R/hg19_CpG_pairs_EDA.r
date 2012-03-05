# Peter Hickey
# 04/01/2012
# EDA of CpG pairings in hg19

# Load the necesssary packages and functions
source("tabulate_comethylation.r")

# Import the data
hg19_CpGs <- import("/usr/local/work/hickey/CpGs/hg19_CpGs.bed.gz")

# Construct the CpG-pairs
disjoint <- makeCpGPairs(hg19_CpGs, "disjoint")
overlapping <- makeCpGPairs(hg19_CpGs, "overlapping")

# Compute and summarise the intra-pair distances
disjoint.widths <- unlist(lapply(X = disjoint, FUN = width))
overlapping.widths <- unlist(lapply(X = overlapping, FUN = width))
summary(disjoint.widths)
summary(overlapping.widths)

# See how many CpG pairs have an intra-pair distance less the upper bound for investigating comethylation
max.read.length <- 80 # Length of the longest read in the library
sum(disjoint.widths <= (max.read.length - 2))
sum(overlapping.widths <= (max.read.length - 2))

# Subset the CpG-pairs object by excluding those pairs with intra-pair distance > the upper bound
# keepPair() returns TRUE if the CpG-pair is to be retained for the comethlyation analysis
keepPair <- function(x, read.length){
  keepers <- width(x) <= (read.length - 2)
  x[keepers, ]
}
# Set the number of cores to use with mclapply
options("cores" = 10)
disjoint.close <- mclapply(X = disjoint, FUN = keepPair, read.length = max.read.length)
overlapping.close <- mclapply(X = overlapping, FUN = keepPair, read.length = max.read.length)

# Calculate how many CpG-pairs are left (and how many there were to begin with)
sum(unlist(lapply(disjoint, nrow)))
sum(unlist(lapply(disjoint.close, nrow)))
sum(unlist(lapply(overlapping, nrow)))
sum(unlist(lapply(overlapping.close, nrow)))