# Compute and plot the SVD of the log odds ratios of comethylation for the Illumina data from Hansen et al.
# Peter Hickey
# 21/12/2011

# Load the data
load("~/Desktop/Comethylation/R/lor/SRR206931_lor.RData")
load("~/Desktop/Comethylation/R/lor/SRR206932_lor.RData")
load("~/Desktop/Comethylation/R/lor/SRR206933_lor.RData")
load("~/Desktop/Comethylation/R/lor/SRR206934_lor.RData")
load("~/Desktop/Comethylation/R/lor/SRR206935_lor.RData")
load("~/Desktop/Comethylation/R/lor/SRR207940_lor.RData")

# Load the necessary packages
library(reshape)

# Take a list of log odds ratios (and optionally their standard errors). Each element of the list corresponds to a chromosome.
# Returns the median log odds ratio (LOR) for each lag (aggregation is performed across chromosomes)
medianLOR <- function(x){
  x_melt <- melt(x)
  # Add the distance variable - ideally this should be given in the lor data but this is not yet the case.
  dist <- as.numeric(rownames(x$chr1)) - 1
  x_melt$dist <- rep(dist, 23 * 2)
  # Compute the median log odds ratio for each lag
  median_lor <- cast(data = x_melt, formula = dist ~ ., fun.aggregate = median, subset=variable=="lor")
  # Chop of the last 5 values since the maximum lag is really read_length - 5 = 80 - 4 = 75
  return(median_lor[1:75, ])
}

SRR206931.mlor <- medianLOR(SRR206931.lor)
SRR206932.mlor <- medianLOR(SRR206932.lor)
SRR206933.mlor <- medianLOR(SRR206933.lor)
SRR206934.mlor <- medianLOR(SRR206934.lor)
SRR206935.mlor <- medianLOR(SRR206935.lor)
SRR207940.mlor <- medianLOR(SRR207940.lor)


mlor <- matrix(c(SRR206931.mlor$'(all)', SRR206932.mlor$'(all)', SRR206933.mlor$'(all)', SRR206934.mlor$'(all)', SRR206935.mlor$'(all)',
                 SRR207940.mlor$'(all)'), byrow = TRUE, nrow = 6)

# Crappy line plot of median log odds ratios coloured by sample
plot(x = 0:74, y = mlor[1, ], type = "l", ylim = c(0,6))
lines(x = 0:74, y = mlor[2, ], col = "red")
lines(x = 0:74, y = mlor[3, ], col = "blue")
lines(x = 0:74, y = mlor[4, ], col = "green")
lines(x = 0:74, y = mlor[5, ], col = "orange")
lines(x = 0:74, y = mlor[6, ], col = "grey")
# TODO: A better plot with legends etc.

## Figure out how to plot first two principal components
mlor.pca <- princomp(mlor)
summary(mlor.pca)
plot(mlor.pca)
loadings(mlor.pca)
mlor.pc <- predict(mlor.pca)
eqscplot(mlor.pc[, 1:2], type = "n", xlab = "first principal component", ylab = "second principal component")
text(mlor.pc[, 1:2], rep(1:6, 75))

prcomp(mlor, scale = TRUE)
plot(prcomp(mlor))
summary(prcomp(mlor, scale = TRUE))
biplot(prcomp(mlor, scale = TRUE))
biplot(prcomp(t(mlor), scale = TRUE))