# Peter Hickey
# 10/09/2012
# Plots of within-fragment comethylation
# This is version 2 and is a substantial re-write of the existing code.
# Functions are stored separated in plot_comethylation_v2_functions.r

#### TODOs ####
# Stratify by other genomic features

#### Read in command line arguments ####
args <- commandArgs(TRUE)
# The name of the sample of interest
sample.name <- args[1]
# Path to genomic-features.RData object
path.to.gf <- args[2]
# Path to plot_comethylation_v2_functions.r
path.to.functions <- args[3]

#### Source functions (includes loading several libraries) ####
source(path.to.functions)

#### Load WF files, stored as GRanges objects and as created by annotate_WF_files.r ####
load(paste0(sample.name, '_complete_annotated_WF_objects.RData'))
load(paste0(sample.name, '_annotated_WF_objects.RData'))

#### Load genomic-features.RData, created by running create_genomic-features.r ####
load(file = path.to.gf)

#### Change working directory to plots folder and set a couple of plotting parameters ####
setwd('plots/')
min.cov <- 10 # Minimum coverage required for a CpG-pair to be included in quantile plots.
quantile.ylim <- c(-4, 6) # y-limits for quantile plots
quantile.xlim <- c(0, round(max(width(wf.all.gr)) + 10, -1)) # x-limits for quantile plots. Rounds to uppermost 10, e.g. 341 -> 350

#### Quantile plots of pair-specific log odds ratios ####
tmp <- lorQuantilePlot(wf.all.gr, min.cov = min.cov, sample.name = sample.name, zero.nic = FALSE, pair.choice = 'all', xlim = quantile.xlim, ylim = quantile.ylim)
ggsave(filename = paste0(sample.name, '.wf_quantile_plot.all_pairs.pdf'), width = 16.9, height = 10.5, plot = tmp)
rm(tmp)
tmp <- lorQuantilePlot(wf.zero.nic.all.gr, min.cov = min.cov, sample.name = sample.name, zero.nic = TRUE, pair.choice = 'all', xlim = quantile.xlim, ylim = quantile.ylim)
ggsave(filename = paste0(sample.name, '.wf_quantile_plot.all_pairs.zero_NIC.pdf'), width = 16.9, height = 10.5, plot = tmp)
rm(tmp)
tmp <- lorQuantilePlot(wf.outermost.gr, min.cov = min.cov, sample.name = sample.name, zero.nic = FALSE, pair.choice = 'outermost', xlim = quantile.xlim, ylim = quantile.ylim)
ggsave(filename = paste0(sample.name, '.wf_quantile_plot.outermost_pairs.pdf'), width = 16.9, height = 10.5, plot = tmp)
rm(tmp)
tmp <- lorQuantilePlot(wf.zero.nic.outermost.gr, min.cov = min.cov, sample.name = sample.name, zero.nic = TRUE, pair.choice = 'outermost', xlim = quantile.xlim, ylim = quantile.ylim)
ggsave(filename = paste0(sample.name, '.wf_quantile_plot.outermost_pairs.zero_NIC.pdf'), width = 16.9, height = 10.5, plot = tmp)
rm(tmp)

#### Quantile plots of pair-specific log odds ratios where both CpGs are within a CGI ####
tmp <- lorQuantilePlot(subsetByOverlaps(wf.all.gr, CGI, type = 'within'), min.cov = min.cov, sample.name = sample.name, zero.nic = FALSE, pair.choice = 'all', genomic.feature = 'CpG-pairs within CGI', xlim = quantile.xlim, ylim = quantile.ylim)
ggsave(filename = paste0(sample.name, '.wf_quantile_plot.all_pairs.CpG-pairs_within_CGI.pdf'), width = 16.9, height = 10.5, plot = tmp)
rm(tmp)
tmp <- lorQuantilePlot(subsetByOverlaps(wf.zero.nic.all.gr, CGI, type = 'within'), min.cov = min.cov, sample.name = sample.name, zero.nic = TRUE, pair.choice = 'all', genomic.feature = 'CpG-pairs within CGI', xlim = quantile.xlim, ylim = quantile.ylim)
ggsave(filename = paste0(sample.name, '.wf_quantile_plot.all_pairs.zero_NIC.CpG-pairs_within_CGI.pdf'), width = 16.9, height = 10.5, plot = tmp)
rm(tmp)
tmp <- lorQuantilePlot(subsetByOverlaps(wf.outermost.gr, CGI, type = 'within'), min.cov = min.cov, sample.name = sample.name, zero.nic = FALSE, pair.choice = 'outermost', genomic.feature = 'CpG-pairs within CGI', xlim = quantile.xlim, ylim = quantile.ylim)
ggsave(filename = paste0(sample.name, '.wf_quantile_plot.outermost_pairs.CpG-pairs_within_CGI.pdf'), width = 16.9, height = 10.5, plot = tmp)
rm(tmp)
tmp <- lorQuantilePlot(subsetByOverlaps(wf.zero.nic.outermost.gr, CGI, type = 'within'), min.cov = min.cov, sample.name = sample.name, zero.nic = TRUE, pair.choice = 'outermost', genomic.feature = 'CpG-pairs within CGI', xlim = quantile.xlim, ylim = quantile.ylim)
ggsave(filename = paste0(sample.name, '.wf_quantile_plot.outermost_pairs.zero_NIC.CpG-pairs_within_CGI.pdf'), width = 16.9, height = 10.5, plot = tmp)
rm(tmp)

#### Quantile plots of pair-specific log odds ratios where both CpGs are outside of a CGI ####
tmp <- lorQuantilePlot(subsetByOverlaps(wf.all.gr, gaps(CGI), type = 'within'), min.cov = min.cov, sample.name = sample.name, zero.nic = FALSE, pair.choice = 'all', genomic.feature = 'CpG-pairs outside a CGI', xlim = quantile.xlim, ylim = quantile.ylim)
ggsave(filename = paste0(sample.name, '.wf_quantile_plot.all_pairs.CpG-pairs_outside_a_CGI.pdf'), width = 16.9, height = 10.5, plot = tmp)
rm(tmp)
tmp <- lorQuantilePlot(subsetByOverlaps(wf.zero.nic.all.gr, gaps(CGI), type = 'within'), min.cov = min.cov, sample.name = sample.name, zero.nic = TRUE, pair.choice = 'all', genomic.feature = 'CpG-pairs outside a CGI', xlim = quantile.xlim, ylim = quantile.ylim)
ggsave(filename = paste0(sample.name, '.wf_quantile_plot.all_pairs.zero_NIC.CpG-pairs_outside_a_CGI.pdf'), width = 16.9, height = 10.5, plot = tmp)
rm(tmp)
tmp <- lorQuantilePlot(subsetByOverlaps(wf.outermost.gr, gaps(CGI), type = 'within'), min.cov = min.cov, sample.name = sample.name, zero.nic = FALSE, pair.choice = 'outermost', genomic.feature = 'CpG-pairs outside a CGI', xlim = quantile.xlim, ylim = quantile.ylim)
ggsave(filename = paste0(sample.name, '.wf_quantile_plot.outermost_pairs.CpG-pairs_outside_a_CGI.pdf'), width = 16.9, height = 10.5, plot = tmp)
rm(tmp)
tmp <- lorQuantilePlot(subsetByOverlaps(wf.zero.nic.outermost.gr, gaps(CGI), type = 'within'), min.cov = min.cov, sample.name = sample.name, zero.nic = TRUE, pair.choice = 'outermost', genomic.feature = 'CpG-pairs outside a CGI', xlim = quantile.xlim, ylim = quantile.ylim)
ggsave(filename = paste0(sample.name, '.wf_quantile_plot.outermost_pairs.zero_NIC.CpG-pairs_outside_a_CGI.pdf'), width = 16.9, height = 10.5, plot = tmp)
rm(tmp)

#### Plots of aggregate log odds ratios ####
all.lor <- aggregatedLor(complete.wf.all.gr, ipd = TRUE, nic = FALSE, correct = TRUE)

aggregate.ylim <- c(-4, 6) # Define x-limits
aggregate.xlim <- c(0, round(max(all.lor$IPD) + 10, -1)) # Define y-limits

tmp <- aggregateLorPlot(all.lor, sample.name, zero.nic = FALSE, pair.choice = 'all', xlim = aggregate.xlim, ylim = aggregate.ylim)
ggsave(filename = paste0(sample.name, '.wf_aggregated_lor_plot.all_pairs.pdf'), width = 16.9, height = 10.5, plot = tmp)
rm(tmp)
outermost.lor <- aggregatedLor(complete.wf.outermost.gr, ipd = TRUE, nic = FALSE, correct = TRUE)
tmp <- aggregateLorPlot(outermost.lor, sample.name, zero.nic = FALSE, pair.choice = 'outermost', xlim = aggregate.xlim, ylim = aggregate.ylim)
ggsave(filename = paste0(sample.name, '.wf_aggregated_lor_plot.outermost_pairs.pdf'), width = 16.9, height = 10.5, plot = tmp)
rm(tmp)
all.zero.nic.lor <- aggregatedLor(complete.wf.zero.nic.all.gr, ipd = TRUE, nic = FALSE, correct = TRUE)
tmp <- aggregateLorPlot(all.zero.nic.lor, sample.name, zero.nic = TRUE, pair.choice = 'all', xlim = aggregate.xlim, ylim = aggregate.ylim)
ggsave(filename = paste0(sample.name, '.wf_aggregated_lor_plot.all_pairs.zero_NIC.pdf'), width = 16.9, height = 10.5, plot = tmp)
rm(tmp)
outermost.zero.nic.lor <- aggregatedLor(complete.wf.zero.nic.outermost.gr, ipd = TRUE, nic = FALSE, correct = TRUE)
tmp <- aggregateLorPlot(outermost.zero.nic.lor, sample.name, zero.nic = TRUE, pair.choice = 'outermost', xlim = aggregate.xlim, ylim = aggregate.ylim)
ggsave(filename = paste0(sample.name, '.wf_aggregated_lor_plot.outermost_pairs.zero_NIC.pdf'), width = 16.9, height = 10.5, plot = tmp)
rm(tmp)

#### Plots of aggregate log odds ratios where both CpGs are within a CGI ####
all.within.cgi.lor <- aggregatedLor(subsetByOverlaps(complete.wf.all.gr, CGI, type = 'within'), ipd = TRUE, nic = FALSE, correct = TRUE)
tmp <- aggregateLorPlot(all.within.cgi.lor, sample.name, zero.nic = FALSE, pair.choice = 'all', genomic.feature = 'CpG-pairs within CGI', xlim = aggregate.xlim, ylim = aggregate.ylim)
ggsave(filename = paste0(sample.name, '.wf_aggregated_lor_plot.all_pairs.CpG-pairs_within_CGI.pdf'), width = 16.9, height = 10.5, plot = tmp)
rm(tmp)
outermost.within.cgi.lor <- aggregatedLor(subsetByOverlaps(complete.wf.outermost.gr, CGI, type = 'within'), ipd = TRUE, nic = FALSE, correct = TRUE)
tmp <- aggregateLorPlot(outermost.within.cgi.lor, sample.name, zero.nic = FALSE, pair.choice = 'outermost', genomic.feature = 'CpG-pairs within CGI', xlim = aggregate.xlim, ylim = aggregate.ylim)
ggsave(filename = paste0(sample.name, '.wf_aggregated_lor_plot.outermost_pairs.CpG-pairs_within_CGI.pdf'), width = 16.9, height = 10.5, plot = tmp)
rm(tmp)
all.zero.nic.within.cgi.lor <- aggregatedLor(subsetByOverlaps(complete.wf.zero.nic.all.gr, CGI, type = 'within'), ipd = TRUE, nic = FALSE, correct = TRUE)
tmp <- aggregateLorPlot(all.zero.nic.within.cgi.lor, sample.name, zero.nic = TRUE, pair.choice = 'all', genomic.feature = 'CpG-pairs within CGI', xlim = aggregate.xlim, ylim = aggregate.ylim)
ggsave(filename = paste0(sample.name, '.wf_aggregated_lor_plot.all_pairs.zero_NIC.CpG-pairs_within_CGI.pdf'), width = 16.9, height = 10.5, plot = tmp)
rm(tmp)
outermost.zero.nic.within.cgi.lor <- aggregatedLor(subsetByOverlaps(complete.wf.zero.nic.outermost.gr, CGI, type = 'within'), ipd = TRUE, nic = FALSE, correct = TRUE)
tmp <- aggregateLorPlot(outermost.zero.nic.within.cgi.lor, sample.name, zero.nic = TRUE, pair.choice = 'outermost', genomic.feature = 'CpG-pairs within CGI', xlim = aggregate.xlim, ylim = aggregate.ylim)
ggsave(filename = paste0(sample.name, '.wf_aggregated_lor_plot.outermost_pairs.zero_NIC.CpG-pairs_within_CGI.pdf'), width = 16.9, height = 10.5, plot = tmp)
rm(tmp)

#### Plots of aggregate log odds ratios where both CpGs are outside of a CGI ####
all.outside.cgi.lor <- aggregatedLor(subsetByOverlaps(complete.wf.all.gr, gaps(CGI), type = 'within'), ipd = TRUE, nic = FALSE, correct = TRUE)
tmp <- aggregateLorPlot(all.outside.cgi.lor, sample.name, zero.nic = FALSE, pair.choice = 'all', genomic.feature = 'CpG-pairs outside a CGI', xlim = aggregate.xlim, ylim = aggregate.ylim)
ggsave(filename = paste0(sample.name, '.wf_aggregated_lor_plot.all_pairs.CpG-pairs_outside_a_CGI.pdf'), width = 16.9, height = 10.5, plot = tmp)
rm(tmp)
outermost.outside.cgi.lor <- aggregatedLor(subsetByOverlaps(complete.wf.outermost.gr, gaps(CGI), type = 'within'), ipd = TRUE, nic = FALSE, correct = TRUE)
tmp <- aggregateLorPlot(outermost.outside.cgi.lor, sample.name, zero.nic = FALSE, pair.choice = 'outermost', genomic.feature = 'CpG-pairs outside a CGI', xlim = aggregate.xlim, ylim = aggregate.ylim)
ggsave(filename = paste0(sample.name, '.wf_aggregated_lor_plot.outermost_pairs.CpG-pairs_outside_a_CGI.pdf'), width = 16.9, height = 10.5, plot = tmp)
rm(tmp)
all.zero.nic.outside.cgi.lor <- aggregatedLor(subsetByOverlaps(complete.wf.zero.nic.all.gr, gaps(CGI), type = 'within'), ipd = TRUE, nic = FALSE, correct = TRUE)
tmp <- aggregateLorPlot(all.zero.nic.outside.cgi.lor, sample.name, zero.nic = TRUE, pair.choice = 'all', genomic.feature = 'CpG-pairs outside a CGI', xlim = aggregate.xlim, ylim = aggregate.ylim)
ggsave(filename = paste0(sample.name, '.wf_aggregated_lor_plot.all_pairs.zero_NIC.CpG-pairs_outside_a_CGI.pdf'), width = 16.9, height = 10.5, plot = tmp)
rm(tmp)
outermost.zero.nic.outside.cgi.lor <- aggregatedLor(subsetByOverlaps(complete.wf.zero.nic.outermost.gr, gaps(CGI), type = 'within'), ipd = TRUE, nic = FALSE, correct = TRUE)
tmp <- aggregateLorPlot(outermost.zero.nic.outside.cgi.lor, sample.name, zero.nic = TRUE, pair.choice = 'outermost', genomic.feature = 'CpG-pairs outside a CGI', xlim = aggregate.xlim, ylim = aggregate.ylim)
ggsave(filename = paste0(sample.name, '.wf_aggregated_lor_plot.outermost_pairs.zero_NIC.CpG-pairs_outside_a_CGI.pdf'), width = 16.9, height = 10.5, plot = tmp)
rm(tmp)

#### Create data.frame to be used in plots of aggregate log odds ratios stratified by CGI status and average methylation in a 500bp window ####
# Using complete.wf.all.gr
# Add CGI variable with 3 levels (In, Partially.in, Out) to each CpG-pair
fully.in <- countOverlaps(complete.wf.all.gr, CGI, type = 'within')
fully.out <- countOverlaps(complete.wf.all.gr, gaps(CGI), type = 'within')
partially.in <- as.numeric(fully.in == 0 & fully.out == 0)
elementMetadata(complete.wf.all.gr)$CGI <- ifelse(fully.in, 'In', ifelse(fully.out, 'Out', 'Partially.in'))
# Discretise gamma.500 by tertiles
gamma.500.tertiles <- quantile(elementMetadata(complete.wf.all.gr)$gamma.500, c(1/3, 2/3))
elementMetadata(complete.wf.all.gr)$gamma.500.level <- ifelse(elementMetadata(complete.wf.all.gr)$gamma.500 < gamma.500.tertiles[1], 'Low', ifelse(elementMetadata(complete.wf.all.gr)$gamma.500 >= gamma.500.tertiles[2], 'High', 'Intermediate'))
# Compute aggregate log odds ratios for each strata
lor1 <- aggregatedLor(subset(complete.wf.all.gr, elementMetadata(complete.wf.all.gr)$CGI == 'In' & elementMetadata(complete.wf.all.gr)$gamma.500.level == 'Low'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor2 <- aggregatedLor(subset(complete.wf.all.gr, elementMetadata(complete.wf.all.gr)$CGI == 'In' & elementMetadata(complete.wf.all.gr)$gamma.500.level == 'Intermediate'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor3 <- aggregatedLor(subset(complete.wf.all.gr, elementMetadata(complete.wf.all.gr)$CGI == 'In' & elementMetadata(complete.wf.all.gr)$gamma.500.level == 'High'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor4 <- aggregatedLor(subset(complete.wf.all.gr, elementMetadata(complete.wf.all.gr)$CGI == 'Partially.in' & elementMetadata(complete.wf.all.gr)$gamma.500.level == 'Low'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor5 <- aggregatedLor(subset(complete.wf.all.gr, elementMetadata(complete.wf.all.gr)$CGI == 'Partially.in' & elementMetadata(complete.wf.all.gr)$gamma.500.level == 'Intermediate'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor6 <- aggregatedLor(subset(complete.wf.all.gr, elementMetadata(complete.wf.all.gr)$CGI == 'Partially.in' & elementMetadata(complete.wf.all.gr)$gamma.500.level == 'High'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor7 <- aggregatedLor(subset(complete.wf.all.gr, elementMetadata(complete.wf.all.gr)$CGI == 'Out' & elementMetadata(complete.wf.all.gr)$gamma.500.level == 'Low'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor8 <- aggregatedLor(subset(complete.wf.all.gr, elementMetadata(complete.wf.all.gr)$CGI == 'Out' & elementMetadata(complete.wf.all.gr)$gamma.500.level == 'Intermediate'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor9 <- aggregatedLor(subset(complete.wf.all.gr, elementMetadata(complete.wf.all.gr)$CGI == 'Out' & elementMetadata(complete.wf.all.gr)$gamma.500.level == 'High'), ipd = TRUE, nic = FALSE, correct = FALSE)

# Create data.frame of results
complete.wf.all.lor.df <- rbind(lor1, lor2, lor3, lor4, lor5, lor6, lor7, lor8, lor9)
complete.wf.all.lor.df$CGI <- c(rep('Within CGI', nrow(lor1) + nrow(lor2) + nrow(lor3)), rep('Partially in CGI', nrow(lor4) + nrow(lor5) + nrow(lor6)), rep('Outside CGI', nrow(lor7) + nrow(lor8) + nrow(lor9)))
complete.wf.all.lor.df$gamma.500.level <- c(rep('Lowly methylated regions', nrow(lor1)), rep('Intermediately methylated regions', nrow(lor2)), rep('Highly methylated regions', nrow(lor3)), rep('Lowly methylated regions', nrow(lor4)), rep('Intermediately methylated regions', nrow(lor5)), rep('Highly methylated regions', nrow(lor6)), rep('Lowly methylated regions', nrow(lor7)), rep('Intermediately methylated regions', nrow(lor8)), rep('Highly methylated regions', nrow(lor9)))

# Using complete.wf.outermost.gr #
# Add CGI variable with 3 levels (In, Partioutermosty.in, Out) to each CpG-pair
fully.in <- countOverlaps(complete.wf.outermost.gr, CGI, type = 'within')
fully.out <- countOverlaps(complete.wf.outermost.gr, gaps(CGI), type = 'within')
partioutermosty.in <- as.numeric(fully.in == 0 & fully.out == 0)
elementMetadata(complete.wf.outermost.gr)$CGI <- ifelse(fully.in, 'In', ifelse(fully.out, 'Out', 'Partioutermosty.in'))
# Discretise gamma.500 by tertiles
gamma.500.tertiles <- quantile(elementMetadata(complete.wf.outermost.gr)$gamma.500, c(1/3, 2/3))
elementMetadata(complete.wf.outermost.gr)$gamma.500.level <- ifelse(elementMetadata(complete.wf.outermost.gr)$gamma.500 < gamma.500.tertiles[1], 'Low', ifelse(elementMetadata(complete.wf.outermost.gr)$gamma.500 >= gamma.500.tertiles[2], 'High', 'Intermediate'))
# Compute aggregate log odds ratios for each strata
lor1 <- aggregatedLor(subset(complete.wf.outermost.gr, elementMetadata(complete.wf.outermost.gr)$CGI == 'In' & elementMetadata(complete.wf.outermost.gr)$gamma.500.level == 'Low'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor2 <- aggregatedLor(subset(complete.wf.outermost.gr, elementMetadata(complete.wf.outermost.gr)$CGI == 'In' & elementMetadata(complete.wf.outermost.gr)$gamma.500.level == 'Intermediate'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor3 <- aggregatedLor(subset(complete.wf.outermost.gr, elementMetadata(complete.wf.outermost.gr)$CGI == 'In' & elementMetadata(complete.wf.outermost.gr)$gamma.500.level == 'High'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor4 <- aggregatedLor(subset(complete.wf.outermost.gr, elementMetadata(complete.wf.outermost.gr)$CGI == 'Partioutermosty.in' & elementMetadata(complete.wf.outermost.gr)$gamma.500.level == 'Low'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor5 <- aggregatedLor(subset(complete.wf.outermost.gr, elementMetadata(complete.wf.outermost.gr)$CGI == 'Partioutermosty.in' & elementMetadata(complete.wf.outermost.gr)$gamma.500.level == 'Intermediate'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor6 <- aggregatedLor(subset(complete.wf.outermost.gr, elementMetadata(complete.wf.outermost.gr)$CGI == 'Partioutermosty.in' & elementMetadata(complete.wf.outermost.gr)$gamma.500.level == 'High'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor7 <- aggregatedLor(subset(complete.wf.outermost.gr, elementMetadata(complete.wf.outermost.gr)$CGI == 'Out' & elementMetadata(complete.wf.outermost.gr)$gamma.500.level == 'Low'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor8 <- aggregatedLor(subset(complete.wf.outermost.gr, elementMetadata(complete.wf.outermost.gr)$CGI == 'Out' & elementMetadata(complete.wf.outermost.gr)$gamma.500.level == 'Intermediate'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor9 <- aggregatedLor(subset(complete.wf.outermost.gr, elementMetadata(complete.wf.outermost.gr)$CGI == 'Out' & elementMetadata(complete.wf.outermost.gr)$gamma.500.level == 'High'), ipd = TRUE, nic = FALSE, correct = FALSE)

# Using complete.wf.zero.nic.all.gr #
# Add CGI variable with 3 levels (In, Partioutermosty.in, Out) to each CpG-pair
fully.in <- countOverlaps(complete.wf.zero.nic.all.gr, CGI, type = 'within')
fully.out <- countOverlaps(complete.wf.zero.nic.all.gr, gaps(CGI), type = 'within')
partizero.nic.ally.in <- as.numeric(fully.in == 0 & fully.out == 0)
elementMetadata(complete.wf.zero.nic.all.gr)$CGI <- ifelse(fully.in, 'In', ifelse(fully.out, 'Out', 'Partizero.nic.ally.in'))
# Discretise gamma.500 by tertiles
gamma.500.tertiles <- quantile(elementMetadata(complete.wf.zero.nic.all.gr)$gamma.500, c(1/3, 2/3))
elementMetadata(complete.wf.zero.nic.all.gr)$gamma.500.level <- ifelse(elementMetadata(complete.wf.zero.nic.all.gr)$gamma.500 < gamma.500.tertiles[1], 'Low', ifelse(elementMetadata(complete.wf.zero.nic.all.gr)$gamma.500 >= gamma.500.tertiles[2], 'High', 'Intermediate'))
# Compute aggregate log odds ratios for each strata
lor1 <- aggregatedLor(subset(complete.wf.zero.nic.all.gr, elementMetadata(complete.wf.zero.nic.all.gr)$CGI == 'In' & elementMetadata(complete.wf.zero.nic.all.gr)$gamma.500.level == 'Low'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor2 <- aggregatedLor(subset(complete.wf.zero.nic.all.gr, elementMetadata(complete.wf.zero.nic.all.gr)$CGI == 'In' & elementMetadata(complete.wf.zero.nic.all.gr)$gamma.500.level == 'Intermediate'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor3 <- aggregatedLor(subset(complete.wf.zero.nic.all.gr, elementMetadata(complete.wf.zero.nic.all.gr)$CGI == 'In' & elementMetadata(complete.wf.zero.nic.all.gr)$gamma.500.level == 'High'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor4 <- aggregatedLor(subset(complete.wf.zero.nic.all.gr, elementMetadata(complete.wf.zero.nic.all.gr)$CGI == 'Partizero.nic.ally.in' & elementMetadata(complete.wf.zero.nic.all.gr)$gamma.500.level == 'Low'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor5 <- aggregatedLor(subset(complete.wf.zero.nic.all.gr, elementMetadata(complete.wf.zero.nic.all.gr)$CGI == 'Partizero.nic.ally.in' & elementMetadata(complete.wf.zero.nic.all.gr)$gamma.500.level == 'Intermediate'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor6 <- aggregatedLor(subset(complete.wf.zero.nic.all.gr, elementMetadata(complete.wf.zero.nic.all.gr)$CGI == 'Partizero.nic.ally.in' & elementMetadata(complete.wf.zero.nic.all.gr)$gamma.500.level == 'High'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor7 <- aggregatedLor(subset(complete.wf.zero.nic.all.gr, elementMetadata(complete.wf.zero.nic.all.gr)$CGI == 'Out' & elementMetadata(complete.wf.zero.nic.all.gr)$gamma.500.level == 'Low'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor8 <- aggregatedLor(subset(complete.wf.zero.nic.all.gr, elementMetadata(complete.wf.zero.nic.all.gr)$CGI == 'Out' & elementMetadata(complete.wf.zero.nic.all.gr)$gamma.500.level == 'Intermediate'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor9 <- aggregatedLor(subset(complete.wf.zero.nic.all.gr, elementMetadata(complete.wf.zero.nic.all.gr)$CGI == 'Out' & elementMetadata(complete.wf.zero.nic.all.gr)$gamma.500.level == 'High'), ipd = TRUE, nic = FALSE, correct = FALSE)

# Create data.frame of results
complete.wf.zero.nic.all.lor.df <- rbind(lor1, lor2, lor3, lor4, lor5, lor6, lor7, lor8, lor9)
complete.wf.zero.nic.all.lor.df$CGI <- c(rep('Within CGI', nrow(lor1) + nrow(lor2) + nrow(lor3)), rep('Partizero.nic.ally in CGI', nrow(lor4) + nrow(lor5) + nrow(lor6)), rep('Outside CGI', nrow(lor7) + nrow(lor8) + nrow(lor9)))
complete.wf.zero.nic.all.lor.df$gamma.500.level <- c(rep('Lowly methylated regions', nrow(lor1)), rep('Intermediately methylated regions', nrow(lor2)), rep('Highly methylated regions', nrow(lor3)), rep('Lowly methylated regions', nrow(lor4)), rep('Intermediately methylated regions', nrow(lor5)), rep('Highly methylated regions', nrow(lor6)), rep('Lowly methylated regions', nrow(lor7)), rep('Intermediately methylated regions', nrow(lor8)), rep('Highly methylated regions', nrow(lor9)))

# Using complete.wf.zero.nic.outermost.gr #
# Add CGI variable with 3 levels (In, Partioutermosty.in, Out) to each CpG-pair
fully.in <- countOverlaps(complete.wf.zero.nic.outermost.gr, CGI, type = 'within')
fully.out <- countOverlaps(complete.wf.zero.nic.outermost.gr, gaps(CGI), type = 'within')
partizero.nic.outermosty.in <- as.numeric(fully.in == 0 & fully.out == 0)
elementMetadata(complete.wf.zero.nic.outermost.gr)$CGI <- ifelse(fully.in, 'In', ifelse(fully.out, 'Out', 'Partizero.nic.outermosty.in'))
# Discretise gamma.500 by tertiles
gamma.500.tertiles <- quantile(elementMetadata(complete.wf.zero.nic.outermost.gr)$gamma.500, c(1/3, 2/3))
elementMetadata(complete.wf.zero.nic.outermost.gr)$gamma.500.level <- ifelse(elementMetadata(complete.wf.zero.nic.outermost.gr)$gamma.500 < gamma.500.tertiles[1], 'Low', ifelse(elementMetadata(complete.wf.zero.nic.outermost.gr)$gamma.500 >= gamma.500.tertiles[2], 'High', 'Intermediate'))
# Compute aggregate log odds ratios for each strata
lor1 <- aggregatedLor(subset(complete.wf.zero.nic.outermost.gr, elementMetadata(complete.wf.zero.nic.outermost.gr)$CGI == 'In' & elementMetadata(complete.wf.zero.nic.outermost.gr)$gamma.500.level == 'Low'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor2 <- aggregatedLor(subset(complete.wf.zero.nic.outermost.gr, elementMetadata(complete.wf.zero.nic.outermost.gr)$CGI == 'In' & elementMetadata(complete.wf.zero.nic.outermost.gr)$gamma.500.level == 'Intermediate'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor3 <- aggregatedLor(subset(complete.wf.zero.nic.outermost.gr, elementMetadata(complete.wf.zero.nic.outermost.gr)$CGI == 'In' & elementMetadata(complete.wf.zero.nic.outermost.gr)$gamma.500.level == 'High'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor4 <- aggregatedLor(subset(complete.wf.zero.nic.outermost.gr, elementMetadata(complete.wf.zero.nic.outermost.gr)$CGI == 'Partizero.nic.outermosty.in' & elementMetadata(complete.wf.zero.nic.outermost.gr)$gamma.500.level == 'Low'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor5 <- aggregatedLor(subset(complete.wf.zero.nic.outermost.gr, elementMetadata(complete.wf.zero.nic.outermost.gr)$CGI == 'Partizero.nic.outermosty.in' & elementMetadata(complete.wf.zero.nic.outermost.gr)$gamma.500.level == 'Intermediate'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor6 <- aggregatedLor(subset(complete.wf.zero.nic.outermost.gr, elementMetadata(complete.wf.zero.nic.outermost.gr)$CGI == 'Partizero.nic.outermosty.in' & elementMetadata(complete.wf.zero.nic.outermost.gr)$gamma.500.level == 'High'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor7 <- aggregatedLor(subset(complete.wf.zero.nic.outermost.gr, elementMetadata(complete.wf.zero.nic.outermost.gr)$CGI == 'Out' & elementMetadata(complete.wf.zero.nic.outermost.gr)$gamma.500.level == 'Low'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor8 <- aggregatedLor(subset(complete.wf.zero.nic.outermost.gr, elementMetadata(complete.wf.zero.nic.outermost.gr)$CGI == 'Out' & elementMetadata(complete.wf.zero.nic.outermost.gr)$gamma.500.level == 'Intermediate'), ipd = TRUE, nic = FALSE, correct = FALSE)
lor9 <- aggregatedLor(subset(complete.wf.zero.nic.outermost.gr, elementMetadata(complete.wf.zero.nic.outermost.gr)$CGI == 'Out' & elementMetadata(complete.wf.zero.nic.outermost.gr)$gamma.500.level == 'High'), ipd = TRUE, nic = FALSE, correct = FALSE)

# Create data.frame of results
complete.wf.zero.nic.outermost.lor.df <- rbind(lor1, lor2, lor3, lor4, lor5, lor6, lor7, lor8, lor9)
complete.wf.zero.nic.outermost.lor.df$CGI <- c(rep('Within CGI', nrow(lor1) + nrow(lor2) + nrow(lor3)), rep('Partizero.nic.outermosty in CGI', nrow(lor4) + nrow(lor5) + nrow(lor6)), rep('Outside CGI', nrow(lor7) + nrow(lor8) + nrow(lor9)))
complete.wf.zero.nic.outermost.lor.df$gamma.500.level <- c(rep('Lowly methylated regions', nrow(lor1)), rep('Intermediately methylated regions', nrow(lor2)), rep('Highly methylated regions', nrow(lor3)), rep('Lowly methylated regions', nrow(lor4)), rep('Intermediately methylated regions', nrow(lor5)), rep('Highly methylated regions', nrow(lor6)), rep('Lowly methylated regions', nrow(lor7)), rep('Intermediately methylated regions', nrow(lor8)), rep('Highly methylated regions', nrow(lor9)))

#### Plots of aggregate log odds ratios stratified by CGI status and average methylation in a 500bp window ####
tmp <- ggplot(data = complete.wf.all.lor.df, aes(x = IPD, y = lor)) + geom_point() + facet_grid(CGI ~ gamma.500.level) + ggtitle(label = paste0(sample.name, ': Within-fragment aggregated comethylation\npair.choice = all')) + scale_x_continuous('Distance between CpGs (bp)') + scale_y_continuous('Log odds ratio') + presentation.theme
ggsave(filename = paste0(sample.name, '.wf_aggregated_lor_plot.all_pairs.3x3.pdf'), width = 16.9, height = 10.5, plot = tmp)
rm(tmp)

tmp <- ggplot(data = complete.wf.outermost.lor.df, aes(x = IPD, y = lor)) + geom_point() + facet_grid(CGI ~ gamma.500.level) + ggtitle(label = paste0(sample.name, ': Within-fragment aggregated comethylation\npair.choice = outermost')) + scale_x_continuous('Distance between CpGs (bp)') + scale_y_continuous('Log odds ratio') + presentation.theme
ggsave(filename = paste0(sample.name, '.wf_aggregated_lor_plot.outermost_pairs.3x3.pdf'), width = 16.9, height = 10.5, plot = tmp)
rm(tmp)

tmp <- ggplot(data = complete.wf.zero.nic.outermost.lor.df, aes(x = IPD, y = lor)) + geom_point() + facet_grid(CGI ~ gamma.500.level) + ggtitle(label = paste0(sample.name, ': Within-fragment aggregated comethylation\npair.choice = zero.nic.outermost, NIC = 0')) + scale_x_continuous('Distance between CpGs (bp)') + scale_y_continuous('Log odds ratio') + presentation.theme
ggsave(filename = paste0(sample.name, '.wf_aggregated_lor_plot.all_pairs.zero_NIC.3x3.pdf'), width = 16.9, height = 10.5, plot = tmp)
rm(tmp)

tmp <- ggplot(data = complete.wf.zero.nic.all.lor.df, aes(x = IPD, y = lor)) + geom_point() + facet_grid(CGI ~ gamma.500.level) + ggtitle(label = paste0(sample.name, ': Within-fragment aggregated comethylation\npair.choice = zero.nic.all, NIC = 0')) + scale_x_continuous('Distance between CpGs (bp)') + scale_y_continuous('Log odds ratio') + presentation.theme
ggsave(filename = paste0(sample.name, '.wf_aggregated_lor_plot.outermost_pairs.zero_NIC.3x3.pdf'), width = 16.9, height = 10.5, plot = tmp)
rm(tmp)

#### Finished ####
setwd('../')
save.image(paste0(sample.name, '_comethylation.RData'))
