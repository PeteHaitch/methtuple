# Peter Hickey
# 06/07/2012
# Aggregated-by-distance within-fragment comethylation plots, stratified by CGI-status and average methylation, plotted by cell-type

#### WARNINGs ####
# Skipping H9_Laurent and HSF1 data due to strong read-position biases in methylation-values

#### TODOs ####
# Check all labels look good
# Save plots to "Powerpoint dimensions"
# Might change axis labels to a darker font colour

#### Load libraries ####
library(ggplot2)

#### Read in the data and construct cell-type-specific data.frames ####
# Somatic
IMR90_r1 <- read.table('~/Lister_2009_BS-seq_data/IMR90/r1/Downloaded_mapped_reads/WF/IMR90_r1_lor_by_methylation_level_aggregated_by_distance.txt', sep = '\t', header = TRUE)
IMR90_r2 <- read.table('~/Lister_2009_BS-seq_data/IMR90/r2/Downloaded_mapped_reads/WF/IMR90_r2_lor_by_methylation_level_aggregated_by_distance.txt', sep = '\t', header = TRUE)
FF <- read.table('~/Lister_2011_BS-seq_data/FF/Downloaded_mapped_reads/WF/FF_lor_by_methylation_level_aggregated_by_distance.txt', sep = '\t', header = TRUE)
ADS <- read.table('~/Lister_2011_BS-seq_data/ADS/Downloaded_mapped_reads/WF/ADS_lor_by_methylation_level_aggregated_by_distance.txt', sep = '\t', header = TRUE)
ADS_adipose <- read.table('~/Lister_2011_BS-seq_data/ADS-adipose/Downloaded_mapped_reads/WF/ADS-adipose_lor_by_methylation_level_aggregated_by_distance.txt', sep = '\t', header = TRUE)

somatic <- rbind(IMR90_r1, IMR90_r2, FF, ADS, ADS_adipose)
somatic$Sample <- c(rep('IMR90_r1', nrow(IMR90_r1)), rep('IMR90_r2', nrow(IMR90_r2)), rep('FF', nrow(FF)), rep('ADS', nrow(ADS)), rep('ADS-adipose', nrow(ADS_adipose)))
somatic$CGI <- ifelse(somatic$CGI == 'CGI', 'CpG island', 'Non CpG island')

# iPSC
ADS_iPSC <- read.table('~/Lister_2011_BS-seq_data/ADS-iPSC/Downloaded_mapped_reads/WF/ADS-iPSC_lor_by_methylation_level_aggregated_by_distance.txt', sep = '\t', header = TRUE)
IMR90_iPSC <- read.table('~/Lister_2011_BS-seq_data/IMR90-iPSC/Downloaded_mapped_reads/WF/IMR90-iPSC_lor_by_methylation_level_aggregated_by_distance.txt', sep = '\t', header = TRUE)
FF_iPSC_6.9 <- read.table('~/Lister_2011_BS-seq_data/FF-iPSC_6.9/Downloaded_mapped_reads/WF/FF-iPSC_6.9_lor_by_methylation_level_aggregated_by_distance.txt', sep = '\t', header = TRUE)
FF_iPSC_19.7 <- read.table('~/Lister_2011_BS-seq_data/FF-iPSC_19.7/Downloaded_mapped_reads/WF/FF-iPSC_19.7_lor_by_methylation_level_aggregated_by_distance.txt', sep = '\t', header = TRUE)
FF_iPSC_19.11 <- read.table('~/Lister_2011_BS-seq_data/FF-iPSC_19.11/Downloaded_mapped_reads/WF/FF-iPSC_19.11_lor_by_methylation_level_aggregated_by_distance.txt', sep = '\t', header = TRUE)

iPSC <- rbind(ADS_iPSC, IMR90_iPSC, FF_iPSC_6.9, FF_iPSC_19.7, FF_iPSC_19.11)
iPSC$Sample <- c(rep('ADS-iPSC', nrow(ADS_iPSC)), rep('IMR90-iPSC', nrow(IMR90_iPSC)), rep('FF-iPSC_6.9', nrow(FF_iPSC_6.9)), rep('FF-iPSC_19.7', nrow(FF_iPSC_19.7)), rep('FF-iPSC_19.11', nrow(FF_iPSC_19.11)))
iPSC$CGI <- ifelse(iPSC$CGI == 'CGI', 'CpG island', 'Non CpG island')

# ES
H1_r1 <- read.table('~/Lister_2009_BS-seq_data/H1/r1/Downloaded_mapped_reads/WF/H1_r1_lor_by_methylation_level_aggregated_by_distance.txt', sep = '\t', header = TRUE)
H1_r2 <- read.table('~/Lister_2009_BS-seq_data/H1/r2/Downloaded_mapped_reads/WF/H1_r2_lor_by_methylation_level_aggregated_by_distance.txt', sep = '\t', header = TRUE)
H9 <- read.table('~/Lister_2011_BS-seq_data/H9/Downloaded_mapped_reads/WF/H9_lor_by_methylation_level_aggregated_by_distance.txt', sep = '\t', header = TRUE)
# Skip H9_Laurent due to strong read-position biases in methylation-values
# H9_Laurent <- read.table('~/Lister_2011_BS-seq_data/H9_Laurent/Downloaded_mapped_reads/WF/H9_Laurent_lor_by_methylation_level_aggregated_by_distance.txt', sep = '\t', header = TRUE)
# Skip HSF1 due to strong read-position biases in methylation-values
# HSF1 <- read.table('~/Lister_2011_BS-seq_data/HSF1/Downloaded_mapped_reads/WF/HSF1_lor_by_methylation_level_aggregated_by_distance.txt', sep = '\t', header = TRUE)

# ES <- rbind(H1_r1, H1_r2, H9, H9_Laurent)
ES <- rbind(H1_r1, H1_r2, H9)
# ES$Sample <- c(rep('H1_r1', nrow(H1_r1)), rep('H1_r2', nrow(H1_r2)), rep('H9', nrow(H9)), rep('H9_Laurent', nrow(H9_Laurent)))
ES$Sample <- c(rep('H1_r1', nrow(H1_r1)), rep('H1_r2', nrow(H1_r2)), rep('H9', nrow(H9)))
ES$CGI <- ifelse(ES$CGI == 'CGI', 'CpG island', 'Non CpG island')

# IVD
H1_plus_BMP4 <- read.table('~/Lister_2011_BS-seq_data/H1+BMP4/Downloaded_mapped_reads/WF/H1+BMP4_lor_by_methylation_level_aggregated_by_distance.txt', sep = '\t', header = TRUE)
FF_iPSC_19.11_plus_BMP4 <- read.table('~/Lister_2011_BS-seq_data/FF-iPSC_19.11+BMP4/Downloaded_mapped_reads/WF/FF-iPSC_19.11+BMP4_lor_by_methylation_level_aggregated_by_distance.txt', sep = '\t', header = TRUE)

IVD <- rbind(H1_plus_BMP4, FF_iPSC_19.11_plus_BMP4)
IVD$Sample <- c(rep('H1+BMP4', nrow(H1_plus_BMP4)), rep('FF-iPSC_19.11+BMP4', nrow(FF_iPSC_19.11_plus_BMP4)))
IVD$CGI <- ifelse(IVD$CGI == 'CGI', 'CpG island', 'Non CpG island')

#### Plots ####
my.opts <- opts(axis.text.x = theme_text(size = 25, colour = 'black'), axis.text.y = theme_text(size = 25, colour = 'black'), axis.title.x = theme_text(size = 30), axis.title.y = theme_text(size = 30, angle = 90, vjust = 0.4), strip.text.x = theme_text(size = 25, face = 'italic'), strip.text.y = theme_text(size = 25, angle = 270, face = 'italic'), plot.title = theme_text(size = 35, face = 'bold'), legend.text = theme_text(size = 20), legend.title = theme_text(size = 25, face = 'bold'))
# Somatic
somatic_plot <- ggplot(data = somatic, aes(x = lag, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma) + opts(title = paste('Somatic cell-types: Aggregated within-fragment comethylation')) 
# Using lines and no transparency. Enforced x- and y-axis limits. Focus on single-end read-lengths.
somatic_plot + geom_line(size = 2) + scale_x_continuous('Distance between CpGs (bp)', limits = c(0, 90)) + scale_y_continuous('Log odds ratio', limits = c(-1, 4.5)) + my.opts
# Using lines and no transparency. Enforced x- and y-axis limits. Focus on (reliable) paired-end read-lengths.
somatic_plot + geom_line(size = 2) + scale_x_continuous('Distance between CpGs (bp)', limits = c(0, 230)) + scale_y_continuous('Log odds ratio', limits = c(-1, 4.5)) + my.opts

# iPSC
iPSC_plot <- ggplot(data = iPSC, aes(x = lag, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma) + opts(title = paste('iPSC cell-types: Within-fragment aggregated comethylation')) 
# Using lines and no transparency. Enforced x- and y-axis limits. Focus on single-end read-lengths.
iPSC_plot + geom_line(size = 2) + scale_x_continuous('Distance between CpGs (bp)', limits = c(0, 90)) + scale_y_continuous('Log odds ratio', limits = c(-1, 4.5)) + my.opts
# Using lines and no transparency. Enforced x- and y-axis limits. Focus on (reliable) paired-end read-lengths.
iPSC_plot + geom_line(size = 2) + scale_x_continuous('Distance between CpGs (bp)', limits = c(0, 230)) + scale_y_continuous('Log odds ratio', limits = c(-1, 4.5)) + my.opts

# ES
ES_plot <- ggplot(data = ES, aes(x = lag, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma) + opts(title = paste('ES-cell cell-types: Aggregated within-fragment comethylation')) 
# Using lines and no transparency. Enforced x- and y-axis limits. Focus on single-end read-lengths.
ES_plot + geom_line(size = 2) + scale_x_continuous('Distance between CpGs (bp)', limits = c(0, 90)) + scale_y_continuous('Log odds ratio', limits = c(-1, 4.5)) + my.opts
# Using lines and no transparency. Enforced x- and y-axis limits. Focus on (reliable) paired-end read-lengths.
ES_plot + geom_line(size = 2) + scale_x_continuous('Distance between CpGs (bp)', limits = c(0, 230)) + scale_y_continuous('Log odds ratio', limits = c(-1, 4.5)) + my.opts

# IVD
IVD_plot <- ggplot(data = IVD, aes(x = lag, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma) + opts(title = paste('IVD cell-types: Within-fragment aggregated comethylation')) 
# Using lines and no transparency. Enforced x- and y-axis limits. Focus on single-end read-lengths.
IVD_plot + geom_line(size = 2) + scale_x_continuous('Distance between CpGs (bp)', limits = c(0, 90)) + scale_y_continuous('Log odds ratio', limits = c(-1, 4.5)) + my.opts
# Using lines and no transparency. Enforced x- and y-axis limits. Focus on (reliable) paired-end read-lengths.
IVD_plot + geom_line(size = 2) + scale_x_continuous('Distance between CpGs (bp)', limits = c(0, 230)) + scale_y_continuous('Log odds ratio', limits = c(-1, 4.5)) + my.opts

#### Finished ####
