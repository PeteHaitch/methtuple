# Compute odds-ratios for lag 1,...,readLength based on output of SAM2MS.py
# Peter Hickey
# 5/3/2012

# Function to compute the log odds-ratio for a comethylation table
## TODO: Check format of x before proceeding
logOdds <- function(x){
  aa <- x["Z", "Z"] + 0.5
  bb <- x["Z", "z"] + 0.5
  cc <- x["z", "Z"] + 0.5
  dd <- x["z", "z"] + 0.5
  lor <- log((aa*dd)/(bb*cc))
  ASE <- sqrt(1/aa + 1/bb + 1/cc + 1/dd) # See pp54 of Agresti
  return(c(lor, ASE))
}

x <- read.table("SRR094471.out", comment.char = "", stringsAsFactors = F, header = F, sep = " ")
colnames(x) <- c("chr", "pos1", "pos2", "X1", "X2")
x$dist <- x$pos2 - x$pos1

# As a function and parallelised
foo <- function(lag, data){
  xx <- subset(data, dist == lag)
  tt <- table(xx[, 4:5])
  logOdds(tt)
}
library(multicore)
lors <- mclapply(X = 2:max(x$dist), FUN = foo, data = x, mc.cores = 10)

library(ggplot2)
df <- data.frame(distance = 2:max(x$dist))
y <- unlist(lors)
df$lor <- y[seq(1, length(y), 2)]
df$ASE <- y[seq(2, length(y), 2)]

limits <- aes(ymax = lor + 2*ASE, ymin = lor - 2*ASE)
m <- ggplot(data = df, aes(x = distance, y = lor)) + geom_point() + geom_errorbar(limits) + opts(title = paste("FF: SRR094471\n", "n = ", nrow(x))) + scale_y_continuous(limits = c(0, 5))
ggsave(filename = "SRR094471_lors.pdf", plot = m)
dev.off()

# Compare to CGIs definied by Rafa et al.
cgi <- read.table("/usr/local/work/hickey/CpGs/Rafas-cpg-islands-hg19.txt", header = T, stringsAsFactors = F)
library(GenomicRanges)
cgi.gr <- GRanges(seqnames = cgi$chr, IRanges(cgi$start, cgi$end))
x.gr <- GRanges(seqnames = x$chr, IRanges(x$pos1, x$pos2))
olaps <- countOverlaps(x.gr, cgi.gr)
save.image(file = "lor_v2.RData")

x$CGI <- olaps

lors.in.CGI <- mclapply(X = 2:max(x$dist), FUN = foo, data = subset(x, CGI==1), mc.cores = 10)
df.in.CGI <- data.frame(distance = 2:max(x$dist))
y <- unlist(lors.in.CGI)
df.in.CGI$lor <- y[seq(1, length(y), 2)]
df.in.CGI$ASE <- y[seq(2, length(y), 2)]
limits <- aes(ymax = lor + 2*ASE, ymin = lor - 2*ASE)
m.in.CGI <- ggplot(data = df.in.CGI, aes(x = distance, y = lor)) + geom_point() + geom_errorbar(limits) + opts(title = paste("FF: SRR094471 in CGI\n", "n = ", sum(x$CGI))) + scale_y_continuous(limits = c(0, 5))
ggsave(filename = "SRR094471_lors_in_CGI.pdf", plot = m.in.CGI)
dev.off()

lors.out.CGI <- mclapply(X = 2:max(x$dist), FUN = foo, data = subset(x, CGI==0), mc.cores = 10)
df.out.CGI <- data.frame(distance = 2:max(x$dist))
y <- unlist(lors.out.CGI)
df.out.CGI$lor <- y[seq(1, length(y), 2)]
df.out.CGI$ASE <- y[seq(2, length(y), 2)]
limits <- aes(ymax = lor + 2*ASE, ymin = lor - 2*ASE)
m.out.CGI <- ggplot(data = df.out.CGI, aes(x = distance, y = lor)) + geom_point() + geom_errorbar(limits) + opts(title = paste("FF: SRR094471 out CGI\n", "n = ", nrow(x) - sum(x$CGI))) + scale_y_continuous(limits = c(0, 5))
ggsave(filename = "SRR094471_lors_out_CGI.pdf", plot = m.out.CGI)
dev.off()
