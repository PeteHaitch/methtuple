###################################################################################################
# Explorations of distance between neighbouring CpGs
###################################################################################################
## Create GRanges object of CpGs in hg19
library("BSgenome.Hsapiens.UCSC.hg19")
params <- new("BSParams", X = Hsapiens, FUN = matchPattern, simplify = TRUE)
cpgs <- bsapply(params, pattern = "CG")

# Function to convert cpgs (an XStringViews object) to a GRanges object
cpgsdf2GR <- function(x){
  GRanges(seqnames = x, ranges = IRanges(start = cpgs[[x]]@start, width = cpgs[[x]]@width))
}
cpgs.GR <- lapply(X = names(cpgs), FUN = cpgsdf2GR)
names(cpgs.GR) <- names(cpgs)
rm(cpgs)

## Calcualte the distance between neighbouring CpGs
# Actually reports the number of bases between Cs in CpGs
cpgDistance <- function(x){
  n <- length(cpgs.GR[[x]])
  distance(cpgs.GR[[x]][1:(n-1)]@ranges, cpgs.GR[[x]][2:n]@ranges) + 1 
}

distances <- lapply(X = names(cpgs.GR), FUN = cpgDistance)
d <- unlist(distances)

## Draw histograms of distances between neighbouring CpGs - both per chromosome and pooled
library(ggplot2)
d.1000 <- d[d < 1000] # d.1000 captures > 99% of CpGs
length(d.1000) / length(d) * 100
# qplot(d.1000, geom = "histogram", binwidth = 100) # Painfully slow. DON'T RUN

#ggplot(data.frame(d=dd),aes(x=d))+
 # geom_histogram(colour="gray",binwidth=5,
  #               position="identity")+theme_bw()

## Draw density plots of distance
plot(density(d.1000))
plot(density(d.1000), xlim = c(0, 100))

# A hack to get a quick histogram
d.1000.table <- table(d.1000)
plot(d.1000.table/sum(d.1000.table), type = "h", main = "Hack histogram of CpG distances (within 1000bp)",
     xlab = "Distance", ylab = "Pseudo-proportion", col = "orange")
lines(density(d.1000), col = "dodgerBlue", lwd = 2)

## Select two arbitrary proximal CpGs. 
# Selection done by eyeballing data in IGV to find candidate region
my.region <- GRanges(seqnames="chr7", ranges=IRanges(start=112091564, end = 112091713))
my.cpgs <- intersect(cpgs.GR$chr7@ranges, my.region@ranges)
# This region contains two CpGs, 5bp apart with high coverage - but almost all sites are methylated

# chr12:132,701,314-132,701,439
my.region <- GRanges(seqnames="chr12", ranges=IRanges(start=132701314, end = 132701439))
intersect(cpgs.GR$chr7@ranges, my.region@ranges)