# Peter Hickey
# 04/04/2012
# Comparison of Pearson's correlation, r^2 (as used in Shoemaker et al. 2010) and log-odds ratio

# Pearson's correlation of the binary vectors
# cor(x, y)

# Function to compute the r^2 value from Shoemaker et al. 2010. Analogous to LD
## This is just (cor(x, y))^2, hence r^2 !!!
r2 <- function(x, y){
  FMM <- sum(x == 1 & y == 1)
  FUU <- sum(x == 0 & y == 0)
  FMU <- sum(x == 1 & y == 0)
  FUM <- sum(x == 0 & y == 1)
  FMstar <- sum(x)
  FUstar <- sum(!x)
  FstarU <- sum(!y)
  FstarM <- sum(y)
  r2 <- (FMM*FUU - FMU*FUM)^2/(FMstar*FUstar*FstarU*FstarM)
  return(r2)
}

# Log-odds ratio
lor <- function(x, y){
  z <- table(x, y)
  lor <- log((z[1,1]*z[2,2])/(z[1,2]*z[2,1]))
  return(lor)
}

# Compare the three measures on x and y
comp <- function(x, y){
  data.frame(cor = cor(x, y), r2 = r2(x, y), lor = lor(x, y))
}

# Data are the from the two leftmost CpGs in Figure 1a from Shoemaker et al. 2010
# x1 is the leftmost CpG
# y1 is the second-leftmost CpG
# 1 = methylated, 0 = unmethylated
x1 <- c(1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1)
y1 <- c(0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1)

## Compute the three measures on the data from Shoemaker et al. 2010
comp(x1, y1)

## Q. What happens if the order of the vectors is reverse?
## A. Results are identical
x2 <- y1
y2 <- x1
comp(x2, y2)

## Q. What happens if the same vector is repeated 10 times - i.e. increase the counts by a factor of 10?
## A. All three measures are stable with respect to scalar multiplication of the counts
x3 <- rep(x1, 10)
y3 <- rep(y1, 10)
comp(x3, y3)

## Q. What happens if the "opposite" methylation states are reported?
## A. All threa measures are stable with respect to opposite states
x4 <- !x1
y4 <- !y1
comp(x4, y4)

## Q. Compute the three measures for two perfectly dependent vectors
## A. All give the strongest possible measure of association
x5 <- x1
y5 <- x5
comp(x5, y5)

## Q. Compute the three measures for two vectors of IID RVs where P(X = 1) = P(X = 0) = 0.5
## A. All show no association
x6 <- rnorm(200, 0, 1) < 0
y6 <- rnorm(200, 0, 1) < 0 
comp(x6, y6)

## Q. Compute the three measures for two vectors of IID RVs where P(X = 1) = P(X = 0) = 0.5 AND the Pearson correlation is negative.
## Both the cor and lor measures are negative giving an indication of the directino of the association whereas the r^2 gives no information on the directionality.
set.seed(100)
x7 <- rnorm(100, 0, 1) < 0
y7 <- rnorm(100, 0, 1) < 0 
comp(x7, y7)

## Q. Compute the three measures for two "opposite" vectors - still entirely dependent/
## A. All give the strongest possible measure of association
x8 <- x1
y8 <- !x1
comp(x8, y8)