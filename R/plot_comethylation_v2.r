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

#### Finished ####
