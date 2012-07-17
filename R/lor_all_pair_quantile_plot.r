library(ggplot2)
library(reshape2)
library(GenomicRanges)

lor.df <- data.frame(lag = width(WF.gr) - 1, lor = elementMetadata(WF.gr)$lor, cov = elementMetadata(WF.gr)$cov)
a <- data.frame(lag = 2:349, Q10 = NA, Q25 = NA, Q50 = NA, Q75 = NA, Q90 = NA )
for(i in 2:349){
  tmp <- subset(lor.df, lor.df$lag == i)
  a[i-1, ] <- c(i, quantile(tmp$lor, c(0.1, 0.25, 0.5, 0.75, 0.9)))
}
b <- melt(a, "lag")

qplot(x = lag, y = value, lty = variable, data = b, xlim = c(0, 240), geom = "line", main = "ADS: Within-fragment comethylation", ylim = c(-2, 5), ylab = "log-odds ratio", xlab = "Distance between CpGs (bp)") + scale_linetype_manual(name="Quantile",  values=c(Q10 = 3, Q25 = 2, Q50 = 1, Q75 = 2, Q90 = 3))
ggsave(filename = 'ADS_all_pairs_lor_quantiles.pdf')

names(lor.df) <- c('lag', 'lor', 'cov')
lor.df <- subset(lor.df, lor.df$cov >= 10)
a <- data.frame(lag = 2:349, Q10 = NA, Q25 = NA, Q50 = NA, Q75 = NA, Q90 = NA)
for(i in 2:349){
  tmp <- subset(lor.df, lor.df$lag == i)
  a[i-1, ] <- c(i, quantile(tmp$lor, c(0.1, 0.25, 0.5, 0.75, 0.9)))
}
b <- melt(a, "lag")

qplot(x = lag, y = value, lty = variable, data = b, geom = "line", xlim = c(0, 240), main = "ADS: Within-fragment comethylation - min cov = 10", ylim = c(-2, 5), ylab = "log-odds ratio", xlab = "Distance between CpGs (bp)") + scale_linetype_manual(name="Quantile",  values=c(Q10 = 3, Q25 = 2, Q50 = 1, Q75 = 2, Q90 = 3))
ggsave(filename = 'ADS_all_pairs_lor_quantiles_mincov=10.pdf')

lor.df <- subset(lor.df, lor.df$cov >= 20)
a <- data.frame(lag = 2:349, Q10 = NA, Q25 = NA, Q50 = NA, Q75 = NA, Q90 = NA)
for(i in 2:349){
  tmp <- subset(lor.df, lor.df$lag == i)
  a[i-1, ] <- c(i, quantile(tmp$lor, c(0.1, 0.25, 0.5, 0.75, 0.9)))
}
b <- melt(a, "lag")

qplot(x = lag, y = value, lty = variable, data = b, geom = "line", xlim = c(0, 240), main = "ADS: Within-fragment comethylation - min cov = 20", ylim = c(-2, 5), ylab = "log-odds ratio", xlab = "Distance between CpGs (bp)") + scale_linetype_manual(name="Quantile",  values=c(Q10 = 3, Q25 = 2, Q50 = 1, Q75 = 2, Q90 = 3))
ggsave(filename = 'ADS_all_pairs_lor_quantiles_mincov=20.pdf')