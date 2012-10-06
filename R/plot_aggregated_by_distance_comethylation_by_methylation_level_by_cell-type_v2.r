# Peter Hickey
# 02/10/2012
# Aggregated-by-distance within-fragment comethylation plotted by cell-type

#### TODOs ####

#### Load libraries ####
library(ggplot2)
library(stringr)

#### Source plot_comethylation_v2_functions.r so that I can use 'presentation.opts' ####
source('~/Comethylation_scripts/plot_comethylation_v2_functions.r')

#### Load all RData files and extract relevant objects ####
H1_r1.env <- new.env()
load('Lister_2009_BS-seq_data/H1/r1/Downloaded_mapped_reads/WF/H1_r1_comethylation.RData', envir = H1_r1.env)
H1_r1.all.lor <- H1_r1.env$all.lor
H1_r1.all.outside.cgi.lor <- H1_r1.env$all.outside.cgi.lor
H1_r1.all.within.cgi.lor <- H1_r1.env$all.within.cgi.lor
H1_r1.all.zero.nic.lor <- H1_r1.env$all.zero.nic.lor
H1_r1.all.zero.nic.outside.cgi.lor <- H1_r1.env$all.zero.nic.outside.cgi.lor
H1_r1.all.zero.nic.within.cgi.lor <- H1_r1.env$all.zero.nic.within.cgi.lor
H1_r1.outermost.lor <- H1_r1.env$outermost.lor
H1_r1.outermost.outside.cgi.lor <- H1_r1.env$outermost.outside.cgi.lor
H1_r1.outermost.within.cgi.lor <- H1_r1.env$outermost.within.cgi.lor
H1_r1.outermost.zero.nic.lor <- H1_r1.env$outermost.zero.nic.lor
H1_r1.outermost.zero.nic.outside.cgi.lor <- H1_r1.env$outermost.zero.nic.outside.cgi.lor
H1_r1.outermost.zero.nic.within.cgi.lor <- H1_r1.env$outermost.zero.nic.within.cgi.lor
H1_r1.complete.wf.all.lor.df <- H1_r1.env$complete.wf.all.lor.df
H1_r1.complete.wf.outermost.lor.df <- H1_r1.env$complete.wf.outermost.lor.df
H1_r1.complete.wf.zero.nic.all.lor.df <- H1_r1.env$complete.wf.zero.nic.all.lor.df
H1_r1.complete.wf.zero.nic.outermost.lor.df <- H1_r1.env$complete.wf.zero.nic.outermost.lor.df
rm(H1_r1.env)
gc()

H1_r2.env <- new.env()
load('Lister_2009_BS-seq_data/H1/r2/Downloaded_mapped_reads/WF/H1_r2_comethylation.RData', envir = H1_r2.env)
H1_r2.all.lor <- H1_r2.env$all.lor
H1_r2.all.outside.cgi.lor <- H1_r2.env$all.outside.cgi.lor
H1_r2.all.within.cgi.lor <- H1_r2.env$all.within.cgi.lor
H1_r2.all.zero.nic.lor <- H1_r2.env$all.zero.nic.lor
H1_r2.all.zero.nic.outside.cgi.lor <- H1_r2.env$all.zero.nic.outside.cgi.lor
H1_r2.all.zero.nic.within.cgi.lor <- H1_r2.env$all.zero.nic.within.cgi.lor
H1_r2.outermost.lor <- H1_r2.env$outermost.lor
H1_r2.outermost.outside.cgi.lor <- H1_r2.env$outermost.outside.cgi.lor
H1_r2.outermost.within.cgi.lor <- H1_r2.env$outermost.within.cgi.lor
H1_r2.outermost.zero.nic.lor <- H1_r2.env$outermost.zero.nic.lor
H1_r2.outermost.zero.nic.outside.cgi.lor <- H1_r2.env$outermost.zero.nic.outside.cgi.lor
H1_r2.outermost.zero.nic.within.cgi.lor <- H1_r2.env$outermost.zero.nic.within.cgi.lor
H1_r2.complete.wf.all.lor.df <- H1_r2.env$complete.wf.all.lor.df
H1_r2.complete.wf.outermost.lor.df <- H1_r2.env$complete.wf.outermost.lor.df
H1_r2.complete.wf.zero.nic.all.lor.df <- H1_r2.env$complete.wf.zero.nic.all.lor.df
H1_r2.complete.wf.zero.nic.outermost.lor.df <- H1_r2.env$complete.wf.zero.nic.outermost.lor.df
rm(H1_r2.env)
gc()

H1_merged.env <- new.env()
load('Lister_2009_BS-seq_data/H1/merged/Downloaded_mapped_reads/WF/H1_merged_comethylation.RData', envir = H1_merged.env)
H1_merged.all.lor <- H1_merged.env$all.lor
H1_merged.all.outside.cgi.lor <- H1_merged.env$all.outside.cgi.lor
H1_merged.all.within.cgi.lor <- H1_merged.env$all.within.cgi.lor
H1_merged.all.zero.nic.lor <- H1_merged.env$all.zero.nic.lor
H1_merged.all.zero.nic.outside.cgi.lor <- H1_merged.env$all.zero.nic.outside.cgi.lor
H1_merged.all.zero.nic.within.cgi.lor <- H1_merged.env$all.zero.nic.within.cgi.lor
H1_merged.outermost.lor <- H1_merged.env$outermost.lor
H1_merged.outermost.outside.cgi.lor <- H1_merged.env$outermost.outside.cgi.lor
H1_merged.outermost.within.cgi.lor <- H1_merged.env$outermost.within.cgi.lor
H1_merged.outermost.zero.nic.lor <- H1_merged.env$outermost.zero.nic.lor
H1_merged.outermost.zero.nic.outside.cgi.lor <- H1_merged.env$outermost.zero.nic.outside.cgi.lor
H1_merged.outermost.zero.nic.within.cgi.lor <- H1_merged.env$outermost.zero.nic.within.cgi.lor
H1_merged.complete.wf.all.lor.df <- H1_merged.env$complete.wf.all.lor.df
H1_merged.complete.wf.outermost.lor.df <- H1_merged.env$complete.wf.outermost.lor.df
H1_merged.complete.wf.zero.nic.all.lor.df <- H1_merged.env$complete.wf.zero.nic.all.lor.df
H1_merged.complete.wf.zero.nic.outermost.lor.df <- H1_merged.env$complete.wf.zero.nic.outermost.lor.df
rm(H1_merged.env)
gc()

IMR90_r1.env <- new.env()
load('Lister_2009_BS-seq_data/IMR90/r1/Downloaded_mapped_reads/WF/IMR90_r1_comethylation.RData', envir = IMR90_r1.env)
IMR90_r1.all.lor <- IMR90_r1.env$all.lor
IMR90_r1.all.outside.cgi.lor <- IMR90_r1.env$all.outside.cgi.lor
IMR90_r1.all.within.cgi.lor <- IMR90_r1.env$all.within.cgi.lor
IMR90_r1.all.zero.nic.lor <- IMR90_r1.env$all.zero.nic.lor
IMR90_r1.all.zero.nic.outside.cgi.lor <- IMR90_r1.env$all.zero.nic.outside.cgi.lor
IMR90_r1.all.zero.nic.within.cgi.lor <- IMR90_r1.env$all.zero.nic.within.cgi.lor
IMR90_r1.outermost.lor <- IMR90_r1.env$outermost.lor
IMR90_r1.outermost.outside.cgi.lor <- IMR90_r1.env$outermost.outside.cgi.lor
IMR90_r1.outermost.within.cgi.lor <- IMR90_r1.env$outermost.within.cgi.lor
IMR90_r1.outermost.zero.nic.lor <- IMR90_r1.env$outermost.zero.nic.lor
IMR90_r1.outermost.zero.nic.outside.cgi.lor <- IMR90_r1.env$outermost.zero.nic.outside.cgi.lor
IMR90_r1.outermost.zero.nic.within.cgi.lor <- IMR90_r1.env$outermost.zero.nic.within.cgi.lor
IMR90_r1.complete.wf.all.lor.df <- IMR90_r1.env$complete.wf.all.lor.df
IMR90_r1.complete.wf.outermost.lor.df <- IMR90_r1.env$complete.wf.outermost.lor.df
IMR90_r1.complete.wf.zero.nic.all.lor.df <- IMR90_r1.env$complete.wf.zero.nic.all.lor.df
IMR90_r1.complete.wf.zero.nic.outermost.lor.df <- IMR90_r1.env$complete.wf.zero.nic.outermost.lor.df
rm(IMR90_r1.env)
gc()

IMR90_r2.env <- new.env()
load('Lister_2009_BS-seq_data/IMR90/r2/Downloaded_mapped_reads/WF/IMR90_r2_comethylation.RData', envir = IMR90_r2.env)
IMR90_r2.all.lor <- IMR90_r2.env$all.lor
IMR90_r2.all.outside.cgi.lor <- IMR90_r2.env$all.outside.cgi.lor
IMR90_r2.all.within.cgi.lor <- IMR90_r2.env$all.within.cgi.lor
IMR90_r2.all.zero.nic.lor <- IMR90_r2.env$all.zero.nic.lor
IMR90_r2.all.zero.nic.outside.cgi.lor <- IMR90_r2.env$all.zero.nic.outside.cgi.lor
IMR90_r2.all.zero.nic.within.cgi.lor <- IMR90_r2.env$all.zero.nic.within.cgi.lor
IMR90_r2.outermost.lor <- IMR90_r2.env$outermost.lor
IMR90_r2.outermost.outside.cgi.lor <- IMR90_r2.env$outermost.outside.cgi.lor
IMR90_r2.outermost.within.cgi.lor <- IMR90_r2.env$outermost.within.cgi.lor
IMR90_r2.outermost.zero.nic.lor <- IMR90_r2.env$outermost.zero.nic.lor
IMR90_r2.outermost.zero.nic.outside.cgi.lor <- IMR90_r2.env$outermost.zero.nic.outside.cgi.lor
IMR90_r2.outermost.zero.nic.within.cgi.lor <- IMR90_r2.env$outermost.zero.nic.within.cgi.lor
IMR90_r2.complete.wf.all.lor.df <- IMR90_r2.env$complete.wf.all.lor.df
IMR90_r2.complete.wf.outermost.lor.df <- IMR90_r2.env$complete.wf.outermost.lor.df
IMR90_r2.complete.wf.zero.nic.all.lor.df <- IMR90_r2.env$complete.wf.zero.nic.all.lor.df
IMR90_r2.complete.wf.zero.nic.outermost.lor.df <- IMR90_r2.env$complete.wf.zero.nic.outermost.lor.df
rm(IMR90_r2.env)
gc()

IMR90_merged.env <- new.env()
load('Lister_2009_BS-seq_data/IMR90/merged/Downloaded_mapped_reads/WF/IMR90_merged_comethylation.RData', envir = IMR90_merged.env)
IMR90_merged.all.lor <- IMR90_merged.env$all.lor
IMR90_merged.all.outside.cgi.lor <- IMR90_merged.env$all.outside.cgi.lor
IMR90_merged.all.within.cgi.lor <- IMR90_merged.env$all.within.cgi.lor
IMR90_merged.all.zero.nic.lor <- IMR90_merged.env$all.zero.nic.lor
IMR90_merged.all.zero.nic.outside.cgi.lor <- IMR90_merged.env$all.zero.nic.outside.cgi.lor
IMR90_merged.all.zero.nic.within.cgi.lor <- IMR90_merged.env$all.zero.nic.within.cgi.lor
IMR90_merged.outermost.lor <- IMR90_merged.env$outermost.lor
IMR90_merged.outermost.outside.cgi.lor <- IMR90_merged.env$outermost.outside.cgi.lor
IMR90_merged.outermost.within.cgi.lor <- IMR90_merged.env$outermost.within.cgi.lor
IMR90_merged.outermost.zero.nic.lor <- IMR90_merged.env$outermost.zero.nic.lor
IMR90_merged.outermost.zero.nic.outside.cgi.lor <- IMR90_merged.env$outermost.zero.nic.outside.cgi.lor
IMR90_merged.outermost.zero.nic.within.cgi.lor <- IMR90_merged.env$outermost.zero.nic.within.cgi.lor
IMR90_merged.complete.wf.all.lor.df <- IMR90_merged.env$complete.wf.all.lor.df
IMR90_merged.complete.wf.outermost.lor.df <- IMR90_merged.env$complete.wf.outermost.lor.df
IMR90_merged.complete.wf.zero.nic.all.lor.df <- IMR90_merged.env$complete.wf.zero.nic.all.lor.df
IMR90_merged.complete.wf.zero.nic.outermost.lor.df <- IMR90_merged.env$complete.wf.zero.nic.outermost.lor.df
rm(IMR90_merged.env)
gc()

ADS.env <- new.env()
load('Lister_2011_BS-seq_data/ADS/Downloaded_mapped_reads/WF/ADS_comethylation.RData', envir = ADS.env)
ADS.all.lor <- ADS.env$all.lor
ADS.all.outside.cgi.lor <- ADS.env$all.outside.cgi.lor
ADS.all.within.cgi.lor <- ADS.env$all.within.cgi.lor
ADS.all.zero.nic.lor <- ADS.env$all.zero.nic.lor
ADS.all.zero.nic.outside.cgi.lor <- ADS.env$all.zero.nic.outside.cgi.lor
ADS.all.zero.nic.within.cgi.lor <- ADS.env$all.zero.nic.within.cgi.lor
ADS.outermost.lor <- ADS.env$outermost.lor
ADS.outermost.outside.cgi.lor <- ADS.env$outermost.outside.cgi.lor
ADS.outermost.within.cgi.lor <- ADS.env$outermost.within.cgi.lor
ADS.outermost.zero.nic.lor <- ADS.env$outermost.zero.nic.lor
ADS.outermost.zero.nic.outside.cgi.lor <- ADS.env$outermost.zero.nic.outside.cgi.lor
ADS.outermost.zero.nic.within.cgi.lor <- ADS.env$outermost.zero.nic.within.cgi.lor
ADS.complete.wf.all.lor.df <- ADS.env$complete.wf.all.lor.df
ADS.complete.wf.outermost.lor.df <- ADS.env$complete.wf.outermost.lor.df
ADS.complete.wf.zero.nic.all.lor.df <- ADS.env$complete.wf.zero.nic.all.lor.df
ADS.complete.wf.zero.nic.outermost.lor.df <- ADS.env$complete.wf.zero.nic.outermost.lor.df
rm(ADS.env)
gc()

ADS_adipose.env <- new.env()
load('Lister_2011_BS-seq_data/ADS-adipose/Downloaded_mapped_reads/WF/ADS-adipose_comethylation.RData', envir = ADS_adipose.env)
ADS_adipose.all.lor <- ADS_adipose.env$all.lor
ADS_adipose.all.outside.cgi.lor <- ADS_adipose.env$all.outside.cgi.lor
ADS_adipose.all.within.cgi.lor <- ADS_adipose.env$all.within.cgi.lor
ADS_adipose.all.zero.nic.lor <- ADS_adipose.env$all.zero.nic.lor
ADS_adipose.all.zero.nic.outside.cgi.lor <- ADS_adipose.env$all.zero.nic.outside.cgi.lor
ADS_adipose.all.zero.nic.within.cgi.lor <- ADS_adipose.env$all.zero.nic.within.cgi.lor
ADS_adipose.outermost.lor <- ADS_adipose.env$outermost.lor
ADS_adipose.outermost.outside.cgi.lor <- ADS_adipose.env$outermost.outside.cgi.lor
ADS_adipose.outermost.within.cgi.lor <- ADS_adipose.env$outermost.within.cgi.lor
ADS_adipose.outermost.zero.nic.lor <- ADS_adipose.env$outermost.zero.nic.lor
ADS_adipose.outermost.zero.nic.outside.cgi.lor <- ADS_adipose.env$outermost.zero.nic.outside.cgi.lor
ADS_adipose.outermost.zero.nic.within.cgi.lor <- ADS_adipose.env$outermost.zero.nic.within.cgi.lor
ADS_adipose.complete.wf.all.lor.df <- ADS_adipose.env$complete.wf.all.lor.df
ADS_adipose.complete.wf.outermost.lor.df <- ADS_adipose.env$complete.wf.outermost.lor.df
ADS_adipose.complete.wf.zero.nic.all.lor.df <- ADS_adipose.env$complete.wf.zero.nic.all.lor.df
ADS_adipose.complete.wf.zero.nic.outermost.lor.df <- ADS_adipose.env$complete.wf.zero.nic.outermost.lor.df
rm(ADS_adipose.env)
gc()

ADS_iPSC.env <- new.env()
load('Lister_2011_BS-seq_data/ADS_iPSC/Downloaded_mapped_reads/WF/ADS_iPSC_comethylation.RData', envir = ADS_iPSC.env)
ADS_iPSC.all.lor <- ADS_iPSC.env$all.lor
ADS_iPSC.all.outside.cgi.lor <- ADS_iPSC.env$all.outside.cgi.lor
ADS_iPSC.all.within.cgi.lor <- ADS_iPSC.env$all.within.cgi.lor
ADS_iPSC.all.zero.nic.lor <- ADS_iPSC.env$all.zero.nic.lor
ADS_iPSC.all.zero.nic.outside.cgi.lor <- ADS_iPSC.env$all.zero.nic.outside.cgi.lor
ADS_iPSC.all.zero.nic.within.cgi.lor <- ADS_iPSC.env$all.zero.nic.within.cgi.lor
ADS_iPSC.outermost.lor <- ADS_iPSC.env$outermost.lor
ADS_iPSC.outermost.outside.cgi.lor <- ADS_iPSC.env$outermost.outside.cgi.lor
ADS_iPSC.outermost.within.cgi.lor <- ADS_iPSC.env$outermost.within.cgi.lor
ADS_iPSC.outermost.zero.nic.lor <- ADS_iPSC.env$outermost.zero.nic.lor
ADS_iPSC.outermost.zero.nic.outside.cgi.lor <- ADS_iPSC.env$outermost.zero.nic.outside.cgi.lor
ADS_iPSC.outermost.zero.nic.within.cgi.lor <- ADS_iPSC.env$outermost.zero.nic.within.cgi.lor
ADS_iPSC.complete.wf.all.lor.df <- ADS_iPSC.env$complete.wf.all.lor.df
ADS_iPSC.complete.wf.outermost.lor.df <- ADS_iPSC.env$complete.wf.outermost.lor.df
ADS_iPSC.complete.wf.zero.nic.all.lor.df <- ADS_iPSC.env$complete.wf.zero.nic.all.lor.df
ADS_iPSC.complete.wf.zero.nic.outermost.lor.df <- ADS_iPSC.env$complete.wf.zero.nic.outermost.lor.df
rm(ADS_iPSC.env)
gc()

FF.env <- new.env()
load('Lister_2011_BS-seq_data/FF/Downloaded_mapped_reads/WF/FF_comethylation.RData', envir = FF.env)
FF.all.lor <- FF.env$all.lor
FF.all.outside.cgi.lor <- FF.env$all.outside.cgi.lor
FF.all.within.cgi.lor <- FF.env$all.within.cgi.lor
FF.all.zero.nic.lor <- FF.env$all.zero.nic.lor
FF.all.zero.nic.outside.cgi.lor <- FF.env$all.zero.nic.outside.cgi.lor
FF.all.zero.nic.within.cgi.lor <- FF.env$all.zero.nic.within.cgi.lor
FF.outermost.lor <- FF.env$outermost.lor
FF.outermost.outside.cgi.lor <- FF.env$outermost.outside.cgi.lor
FF.outermost.within.cgi.lor <- FF.env$outermost.within.cgi.lor
FF.outermost.zero.nic.lor <- FF.env$outermost.zero.nic.lor
FF.outermost.zero.nic.outside.cgi.lor <- FF.env$outermost.zero.nic.outside.cgi.lor
FF.outermost.zero.nic.within.cgi.lor <- FF.env$outermost.zero.nic.within.cgi.lor
FF.complete.wf.all.lor.df <- FF.env$complete.wf.all.lor.df
FF.complete.wf.outermost.lor.df <- FF.env$complete.wf.outermost.lor.df
FF.complete.wf.zero.nic.all.lor.df <- FF.env$complete.wf.zero.nic.all.lor.df
FF.complete.wf.zero.nic.outermost.lor.df <- FF.env$complete.wf.zero.nic.outermost.lor.df
rm(FF.env)
gc()

FF_iPSC_19.11.env <- new.env()
load('Lister_2011_BS-seq_data/FF-iPSC_19.11/Downloaded_mapped_reads/WF/FF-iPSC_19.11_comethylation.RData', envir = FF_iPSC_19.11.env)
FF_iPSC_19.11.all.lor <- FF_iPSC_19.11.env$all.lor
FF_iPSC_19.11.all.outside.cgi.lor <- FF_iPSC_19.11.env$all.outside.cgi.lor
FF_iPSC_19.11.all.within.cgi.lor <- FF_iPSC_19.11.env$all.within.cgi.lor
FF_iPSC_19.11.all.zero.nic.lor <- FF_iPSC_19.11.env$all.zero.nic.lor
FF_iPSC_19.11.all.zero.nic.outside.cgi.lor <- FF_iPSC_19.11.env$all.zero.nic.outside.cgi.lor
FF_iPSC_19.11.all.zero.nic.within.cgi.lor <- FF_iPSC_19.11.env$all.zero.nic.within.cgi.lor
FF_iPSC_19.11.outermost.lor <- FF_iPSC_19.11.env$outermost.lor
FF_iPSC_19.11.outermost.outside.cgi.lor <- FF_iPSC_19.11.env$outermost.outside.cgi.lor
FF_iPSC_19.11.outermost.within.cgi.lor <- FF_iPSC_19.11.env$outermost.within.cgi.lor
FF_iPSC_19.11.outermost.zero.nic.lor <- FF_iPSC_19.11.env$outermost.zero.nic.lor
FF_iPSC_19.11.outermost.zero.nic.outside.cgi.lor <- FF_iPSC_19.11.env$outermost.zero.nic.outside.cgi.lor
FF_iPSC_19.11.outermost.zero.nic.within.cgi.lor <- FF_iPSC_19.11.env$outermost.zero.nic.within.cgi.lor
FF_iPSC_19.11.complete.wf.all.lor.df <- FF_iPSC_19.11.env$complete.wf.all.lor.df
FF_iPSC_19.11.complete.wf.outermost.lor.df <- FF_iPSC_19.11.env$complete.wf.outermost.lor.df
FF_iPSC_19.11.complete.wf.zero.nic.all.lor.df <- FF_iPSC_19.11.env$complete.wf.zero.nic.all.lor.df
FF_iPSC_19.11.complete.wf.zero.nic.outermost.lor.df <- FF_iPSC_19.11.env$complete.wf.zero.nic.outermost.lor.df
rm(FF_iPSC_19.11.env)
gc()

FF_iPSC_19.11_BMP4.env <- new.env()
load('Lister_2011_BS-seq_data/FF-iPSC_19.11+BMP4/Downloaded_mapped_reads/WF/FF-iPSC_19.11+BMP4_comethylation.RData', envir = FF_iPSC_19.11_BMP4.env)
FF_iPSC_19.11_BMP4.all.lor <- FF_iPSC_19.11_BMP4.env$all.lor
FF_iPSC_19.11_BMP4.all.outside.cgi.lor <- FF_iPSC_19.11_BMP4.env$all.outside.cgi.lor
FF_iPSC_19.11_BMP4.all.within.cgi.lor <- FF_iPSC_19.11_BMP4.env$all.within.cgi.lor
FF_iPSC_19.11_BMP4.all.zero.nic.lor <- FF_iPSC_19.11_BMP4.env$all.zero.nic.lor
FF_iPSC_19.11_BMP4.all.zero.nic.outside.cgi.lor <- FF_iPSC_19.11_BMP4.env$all.zero.nic.outside.cgi.lor
FF_iPSC_19.11_BMP4.all.zero.nic.within.cgi.lor <- FF_iPSC_19.11_BMP4.env$all.zero.nic.within.cgi.lor
FF_iPSC_19.11_BMP4.outermost.lor <- FF_iPSC_19.11_BMP4.env$outermost.lor
FF_iPSC_19.11_BMP4.outermost.outside.cgi.lor <- FF_iPSC_19.11_BMP4.env$outermost.outside.cgi.lor
FF_iPSC_19.11_BMP4.outermost.within.cgi.lor <- FF_iPSC_19.11_BMP4.env$outermost.within.cgi.lor
FF_iPSC_19.11_BMP4.outermost.zero.nic.lor <- FF_iPSC_19.11_BMP4.env$outermost.zero.nic.lor
FF_iPSC_19.11_BMP4.outermost.zero.nic.outside.cgi.lor <- FF_iPSC_19.11_BMP4.env$outermost.zero.nic.outside.cgi.lor
FF_iPSC_19.11_BMP4.outermost.zero.nic.within.cgi.lor <- FF_iPSC_19.11_BMP4.env$outermost.zero.nic.within.cgi.lor
FF_iPSC_19.11_BMP4.complete.wf.all.lor.df <- FF_iPSC_19.11_BMP4.env$complete.wf.all.lor.df
FF_iPSC_19.11_BMP4.complete.wf.outermost.lor.df <- FF_iPSC_19.11_BMP4.env$complete.wf.outermost.lor.df
FF_iPSC_19.11_BMP4.complete.wf.zero.nic.all.lor.df <- FF_iPSC_19.11_BMP4.env$complete.wf.zero.nic.all.lor.df
FF_iPSC_19.11_BMP4.complete.wf.zero.nic.outermost.lor.df <- FF_iPSC_19.11_BMP4.env$complete.wf.zero.nic.outermost.lor.df
rm(FF_iPSC_19.11_BMP4.env)
gc()

FF_iPSC_19.7.env <- new.env()
load('Lister_2011_BS-seq_data/FF-iPSC_19.7/Downloaded_mapped_reads/WF/FF-iPSC_19.7_comethylation.RData', envir = FF_iPSC_19.7.env)
FF_iPSC_19.7.all.lor <- FF_iPSC_19.7.env$all.lor
FF_iPSC_19.7.all.outside.cgi.lor <- FF_iPSC_19.7.env$all.outside.cgi.lor
FF_iPSC_19.7.all.within.cgi.lor <- FF_iPSC_19.7.env$all.within.cgi.lor
FF_iPSC_19.7.all.zero.nic.lor <- FF_iPSC_19.7.env$all.zero.nic.lor
FF_iPSC_19.7.all.zero.nic.outside.cgi.lor <- FF_iPSC_19.7.env$all.zero.nic.outside.cgi.lor
FF_iPSC_19.7.all.zero.nic.within.cgi.lor <- FF_iPSC_19.7.env$all.zero.nic.within.cgi.lor
FF_iPSC_19.7.outermost.lor <- FF_iPSC_19.7.env$outermost.lor
FF_iPSC_19.7.outermost.outside.cgi.lor <- FF_iPSC_19.7.env$outermost.outside.cgi.lor
FF_iPSC_19.7.outermost.within.cgi.lor <- FF_iPSC_19.7.env$outermost.within.cgi.lor
FF_iPSC_19.7.outermost.zero.nic.lor <- FF_iPSC_19.7.env$outermost.zero.nic.lor
FF_iPSC_19.7.outermost.zero.nic.outside.cgi.lor <- FF_iPSC_19.7.env$outermost.zero.nic.outside.cgi.lor
FF_iPSC_19.7.outermost.zero.nic.within.cgi.lor <- FF_iPSC_19.7.env$outermost.zero.nic.within.cgi.lor
FF_iPSC_19.7.complete.wf.all.lor.df <- FF_iPSC_19.7.env$complete.wf.all.lor.df
FF_iPSC_19.7.complete.wf.outermost.lor.df <- FF_iPSC_19.7.env$complete.wf.outermost.lor.df
FF_iPSC_19.7.complete.wf.zero.nic.all.lor.df <- FF_iPSC_19.7.env$complete.wf.zero.nic.all.lor.df
FF_iPSC_19.7.complete.wf.zero.nic.outermost.lor.df <- FF_iPSC_19.7.env$complete.wf.zero.nic.outermost.lor.df
rm(FF_iPSC_19.7.env)
gc()

FF_iPSC_6.9.env <- new.env()
load('Lister_2011_BS-seq_data/FF-iPSC_6.9/Downloaded_mapped_reads/WF/FF-iPSC_6.9_comethylation.RData', envir = FF_iPSC_6.9.env)
FF_iPSC_6.9.all.lor <- FF_iPSC_6.9.env$all.lor
FF_iPSC_6.9.all.outside.cgi.lor <- FF_iPSC_6.9.env$all.outside.cgi.lor
FF_iPSC_6.9.all.within.cgi.lor <- FF_iPSC_6.9.env$all.within.cgi.lor
FF_iPSC_6.9.all.zero.nic.lor <- FF_iPSC_6.9.env$all.zero.nic.lor
FF_iPSC_6.9.all.zero.nic.outside.cgi.lor <- FF_iPSC_6.9.env$all.zero.nic.outside.cgi.lor
FF_iPSC_6.9.all.zero.nic.within.cgi.lor <- FF_iPSC_6.9.env$all.zero.nic.within.cgi.lor
FF_iPSC_6.9.outermost.lor <- FF_iPSC_6.9.env$outermost.lor
FF_iPSC_6.9.outermost.outside.cgi.lor <- FF_iPSC_6.9.env$outermost.outside.cgi.lor
FF_iPSC_6.9.outermost.within.cgi.lor <- FF_iPSC_6.9.env$outermost.within.cgi.lor
FF_iPSC_6.9.outermost.zero.nic.lor <- FF_iPSC_6.9.env$outermost.zero.nic.lor
FF_iPSC_6.9.outermost.zero.nic.outside.cgi.lor <- FF_iPSC_6.9.env$outermost.zero.nic.outside.cgi.lor
FF_iPSC_6.9.outermost.zero.nic.within.cgi.lor <- FF_iPSC_6.9.env$outermost.zero.nic.within.cgi.lor
FF_iPSC_6.9.complete.wf.all.lor.df <- FF_iPSC_6.9.env$complete.wf.all.lor.df
FF_iPSC_6.9.complete.wf.outermost.lor.df <- FF_iPSC_6.9.env$complete.wf.outermost.lor.df
FF_iPSC_6.9.complete.wf.zero.nic.all.lor.df <- FF_iPSC_6.9.env$complete.wf.zero.nic.all.lor.df
FF_iPSC_6.9.complete.wf.zero.nic.outermost.lor.df <- FF_iPSC_6.9.env$complete.wf.zero.nic.outermost.lor.df
rm(FF_iPSC_6.9.env)
gc()

H1_BMP4.env <- new.env()
load('Lister_2011_BS-seq_data/H1+BMP4/Downloaded_mapped_reads/WF/H1+BMP4_comethylation.RData', envir = H1_BMP4.env)
H1_BMP4.all.lor <- H1_BMP4.env$all.lor
H1_BMP4.all.outside.cgi.lor <- H1_BMP4.env$all.outside.cgi.lor
H1_BMP4.all.within.cgi.lor <- H1_BMP4.env$all.within.cgi.lor
H1_BMP4.all.zero.nic.lor <- H1_BMP4.env$all.zero.nic.lor
H1_BMP4.all.zero.nic.outside.cgi.lor <- H1_BMP4.env$all.zero.nic.outside.cgi.lor
H1_BMP4.all.zero.nic.within.cgi.lor <- H1_BMP4.env$all.zero.nic.within.cgi.lor
H1_BMP4.outermost.lor <- H1_BMP4.env$outermost.lor
H1_BMP4.outermost.outside.cgi.lor <- H1_BMP4.env$outermost.outside.cgi.lor
H1_BMP4.outermost.within.cgi.lor <- H1_BMP4.env$outermost.within.cgi.lor
H1_BMP4.outermost.zero.nic.lor <- H1_BMP4.env$outermost.zero.nic.lor
H1_BMP4.outermost.zero.nic.outside.cgi.lor <- H1_BMP4.env$outermost.zero.nic.outside.cgi.lor
H1_BMP4.outermost.zero.nic.within.cgi.lor <- H1_BMP4.env$outermost.zero.nic.within.cgi.lor
H1_BMP4.complete.wf.all.lor.df <- H1_BMP4.env$complete.wf.all.lor.df
H1_BMP4.complete.wf.outermost.lor.df <- H1_BMP4.env$complete.wf.outermost.lor.df
H1_BMP4.complete.wf.zero.nic.all.lor.df <- H1_BMP4.env$complete.wf.zero.nic.all.lor.df
H1_BMP4.complete.wf.zero.nic.outermost.lor.df <- H1_BMP4.env$complete.wf.zero.nic.outermost.lor.df
rm(H1_BMP4.env)
gc()

H9.env <- new.env()
load('Lister_2011_BS-seq_data/H9/Downloaded_mapped_reads/WF/H9_comethylation.RData', envir = H9.env)
H9.all.lor <- H9.env$all.lor
H9.all.outside.cgi.lor <- H9.env$all.outside.cgi.lor
H9.all.within.cgi.lor <- H9.env$all.within.cgi.lor
H9.all.zero.nic.lor <- H9.env$all.zero.nic.lor
H9.all.zero.nic.outside.cgi.lor <- H9.env$all.zero.nic.outside.cgi.lor
H9.all.zero.nic.within.cgi.lor <- H9.env$all.zero.nic.within.cgi.lor
H9.outermost.lor <- H9.env$outermost.lor
H9.outermost.outside.cgi.lor <- H9.env$outermost.outside.cgi.lor
H9.outermost.within.cgi.lor <- H9.env$outermost.within.cgi.lor
H9.outermost.zero.nic.lor <- H9.env$outermost.zero.nic.lor
H9.outermost.zero.nic.outside.cgi.lor <- H9.env$outermost.zero.nic.outside.cgi.lor
H9.outermost.zero.nic.within.cgi.lor <- H9.env$outermost.zero.nic.within.cgi.lor
H9.complete.wf.all.lor.df <- H9.env$complete.wf.all.lor.df
H9.complete.wf.outermost.lor.df <- H9.env$complete.wf.outermost.lor.df
H9.complete.wf.zero.nic.all.lor.df <- H9.env$complete.wf.zero.nic.all.lor.df
H9.complete.wf.zero.nic.outermost.lor.df <- H9.env$complete.wf.zero.nic.outermost.lor.df
rm(H9.env)
gc()

H9_Laurent.env <- new.env()
load('Lister_2011_BS-seq_data/H9_Laurent/Downloaded_mapped_reads/WF/H9_Laurent_comethylation.RData', envir = H9_Laurent.env)
H9_Laurent.all.lor <- H9_Laurent.env$all.lor
H9_Laurent.all.outside.cgi.lor <- H9_Laurent.env$all.outside.cgi.lor
H9_Laurent.all.within.cgi.lor <- H9_Laurent.env$all.within.cgi.lor
H9_Laurent.all.zero.nic.lor <- H9_Laurent.env$all.zero.nic.lor
H9_Laurent.all.zero.nic.outside.cgi.lor <- H9_Laurent.env$all.zero.nic.outside.cgi.lor
H9_Laurent.all.zero.nic.within.cgi.lor <- H9_Laurent.env$all.zero.nic.within.cgi.lor
H9_Laurent.outermost.lor <- H9_Laurent.env$outermost.lor
H9_Laurent.outermost.outside.cgi.lor <- H9_Laurent.env$outermost.outside.cgi.lor
H9_Laurent.outermost.within.cgi.lor <- H9_Laurent.env$outermost.within.cgi.lor
H9_Laurent.outermost.zero.nic.lor <- H9_Laurent.env$outermost.zero.nic.lor
H9_Laurent.outermost.zero.nic.outside.cgi.lor <- H9_Laurent.env$outermost.zero.nic.outside.cgi.lor
H9_Laurent.outermost.zero.nic.within.cgi.lor <- H9_Laurent.env$outermost.zero.nic.within.cgi.lor
H9_Laurent.complete.wf.all.lor.df <- H9_Laurent.env$complete.wf.all.lor.df
H9_Laurent.complete.wf.outermost.lor.df <- H9_Laurent.env$complete.wf.outermost.lor.df
H9_Laurent.complete.wf.zero.nic.all.lor.df <- H9_Laurent.env$complete.wf.zero.nic.all.lor.df
H9_Laurent.complete.wf.zero.nic.outermost.lor.df <- H9_Laurent.env$complete.wf.zero.nic.outermost.lor.df
rm(H9_Laurent.env)
gc()

HSF1.env <- new.env()
load('Lister_2011_BS-seq_data/HSF1/Downloaded_mapped_reads/WF/HSF1_comethylation.RData', envir = HSF1.env)
HSF1.all.lor <- HSF1.env$all.lor
HSF1.all.outside.cgi.lor <- HSF1.env$all.outside.cgi.lor
HSF1.all.within.cgi.lor <- HSF1.env$all.within.cgi.lor
HSF1.all.zero.nic.lor <- HSF1.env$all.zero.nic.lor
HSF1.all.zero.nic.outside.cgi.lor <- HSF1.env$all.zero.nic.outside.cgi.lor
HSF1.all.zero.nic.within.cgi.lor <- HSF1.env$all.zero.nic.within.cgi.lor
HSF1.outermost.lor <- HSF1.env$outermost.lor
HSF1.outermost.outside.cgi.lor <- HSF1.env$outermost.outside.cgi.lor
HSF1.outermost.within.cgi.lor <- HSF1.env$outermost.within.cgi.lor
HSF1.outermost.zero.nic.lor <- HSF1.env$outermost.zero.nic.lor
HSF1.outermost.zero.nic.outside.cgi.lor <- HSF1.env$outermost.zero.nic.outside.cgi.lor
HSF1.outermost.zero.nic.within.cgi.lor <- HSF1.env$outermost.zero.nic.within.cgi.lor
HSF1.complete.wf.all.lor.df <- HSF1.env$complete.wf.all.lor.df
HSF1.complete.wf.outermost.lor.df <- HSF1.env$complete.wf.outermost.lor.df
HSF1.complete.wf.zero.nic.all.lor.df <- HSF1.env$complete.wf.zero.nic.all.lor.df
HSF1.complete.wf.zero.nic.outermost.lor.df <- HSF1.env$complete.wf.zero.nic.outermost.lor.df
rm(HSF1.env)
gc()

IMR90_iPSC.env <- new.env()
load('Lister_2011_BS-seq_data/IMR90-iPSC/Downloaded_mapped_reads/WF/IMR90-iPSC_comethylation.RData', envir = IMR90_iPSC.env)
IMR90_iPSC.all.lor <- IMR90_iPSC.env$all.lor
IMR90_iPSC.all.outside.cgi.lor <- IMR90_iPSC.env$all.outside.cgi.lor
IMR90_iPSC.all.within.cgi.lor <- IMR90_iPSC.env$all.within.cgi.lor
IMR90_iPSC.all.zero.nic.lor <- IMR90_iPSC.env$all.zero.nic.lor
IMR90_iPSC.all.zero.nic.outside.cgi.lor <- IMR90_iPSC.env$all.zero.nic.outside.cgi.lor
IMR90_iPSC.all.zero.nic.within.cgi.lor <- IMR90_iPSC.env$all.zero.nic.within.cgi.lor
IMR90_iPSC.outermost.lor <- IMR90_iPSC.env$outermost.lor
IMR90_iPSC.outermost.outside.cgi.lor <- IMR90_iPSC.env$outermost.outside.cgi.lor
IMR90_iPSC.outermost.within.cgi.lor <- IMR90_iPSC.env$outermost.within.cgi.lor
IMR90_iPSC.outermost.zero.nic.lor <- IMR90_iPSC.env$outermost.zero.nic.lor
IMR90_iPSC.outermost.zero.nic.outside.cgi.lor <- IMR90_iPSC.env$outermost.zero.nic.outside.cgi.lor
IMR90_iPSC.outermost.zero.nic.within.cgi.lor <- IMR90_iPSC.env$outermost.zero.nic.within.cgi.lor
IMR90_iPSC.complete.wf.all.lor.df <- IMR90_iPSC.env$complete.wf.all.lor.df
IMR90_iPSC.complete.wf.outermost.lor.df <- IMR90_iPSC.env$complete.wf.outermost.lor.df
IMR90_iPSC.complete.wf.zero.nic.all.lor.df <- IMR90_iPSC.env$complete.wf.zero.nic.all.lor.df
IMR90_iPSC.complete.wf.zero.nic.outermost.lor.df <- IMR90_iPSC.env$complete.wf.zero.nic.outermost.lor.df
rm(IMR90_iPSC.env)
gc()

save.image('Lister_aggregated_by_distance_comethylation_by_methylation_by_cell-type.RData')

#### Constract somatic data.frames ####
somatic.complete.wf.all.lor.df <- rbind(IMR90_r1.complete.wf.all.lor.df, IMR90_r2.complete.wf.all.lor.df, FF.complete.wf.all.lor.df, ADS.complete.wf.all.lor.df, ADS_adipose.complete.wf.all.lor.df)
somatic.complete.wf.all.lor.df$Sample <- c(rep('IMR90_r1', nrow(IMR90_r1.complete.wf.all.lor.df)), rep('IMR90_r2', nrow(IMR90_r2.complete.wf.all.lor.df)), rep('FF', nrow(FF.complete.wf.all.lor.df)), rep('ADS', nrow(ADS.complete.wf.all.lor.df)), rep('ADS-adipose', nrow(ADS_adipose.complete.wf.all.lor.df)))
# Remove " regions" from gamma.500.level labels to preserve space in plot facet labels
somatic.complete.wf.all.lor.df$gamma.500.level <- str_split_fixed(somatic.complete.wf.all.lor.df$gamma.500.level, ' regions', n = 2)[, 1]

somatic.complete.wf.outermost.lor.df <- rbind(IMR90_r1.complete.wf.outermost.lor.df, IMR90_r2.complete.wf.outermost.lor.df, FF.complete.wf.outermost.lor.df, ADS.complete.wf.outermost.lor.df, ADS_adipose.complete.wf.outermost.lor.df)
somatic.complete.wf.outermost.lor.df$Sample <- c(rep('IMR90_r1', nrow(IMR90_r1.complete.wf.outermost.lor.df)), rep('IMR90_r2', nrow(IMR90_r2.complete.wf.outermost.lor.df)), rep('FF', nrow(FF.complete.wf.outermost.lor.df)), rep('ADS', nrow(ADS.complete.wf.outermost.lor.df)), rep('ADS-adipose', nrow(ADS_adipose.complete.wf.outermost.lor.df)))
# Remove " regions" from gamma.500.level labels to preserve space in plot facet labels
somatic.complete.wf.outermost.lor.df$gamma.500.level <- str_split_fixed(somatic.complete.wf.outermost.lor.df$gamma.500.level, ' regions', n = 2)[, 1]

somatic.complete.wf.zero.nic.all.lor.df <- rbind(IMR90_r1.complete.wf.zero.nic.all.lor.df, IMR90_r2.complete.wf.zero.nic.all.lor.df, FF.complete.wf.zero.nic.all.lor.df, ADS.complete.wf.zero.nic.all.lor.df, ADS_adipose.complete.wf.zero.nic.all.lor.df)
somatic.complete.wf.zero.nic.all.lor.df$Sample <- c(rep('IMR90_r1', nrow(IMR90_r1.complete.wf.zero.nic.all.lor.df)), rep('IMR90_r2', nrow(IMR90_r2.complete.wf.zero.nic.all.lor.df)), rep('FF', nrow(FF.complete.wf.zero.nic.all.lor.df)), rep('ADS', nrow(ADS.complete.wf.zero.nic.all.lor.df)), rep('ADS-adipose', nrow(ADS_adipose.complete.wf.zero.nic.all.lor.df)))
# Remove " regions" from gamma.500.level labels to preserve space in plot facet labels
somatic.complete.wf.zero.nic.all.lor.df$gamma.500.level <- str_split_fixed(somatic.complete.wf.zero.nic.all.lor.df$gamma.500.level, ' regions', n = 2)[, 1]

somatic.complete.wf.zero.nic.outermost.lor.df <- rbind(IMR90_r1.complete.wf.zero.nic.outermost.lor.df, IMR90_r2.complete.wf.zero.nic.outermost.lor.df, FF.complete.wf.zero.nic.outermost.lor.df, ADS.complete.wf.zero.nic.outermost.lor.df, ADS_adipose.complete.wf.zero.nic.outermost.lor.df)
somatic.complete.wf.zero.nic.outermost.lor.df$Sample <- c(rep('IMR90_r1', nrow(IMR90_r1.complete.wf.zero.nic.outermost.lor.df)), rep('IMR90_r2', nrow(IMR90_r2.complete.wf.zero.nic.outermost.lor.df)), rep('FF', nrow(FF.complete.wf.zero.nic.outermost.lor.df)), rep('ADS', nrow(ADS.complete.wf.zero.nic.outermost.lor.df)), rep('ADS-adipose', nrow(ADS_adipose.complete.wf.zero.nic.outermost.lor.df)))
# Remove " regions" from gamma.500.level labels to preserve space in plot facet labels
somatic.complete.wf.zero.nic.outermost.lor.df$gamma.500.level <- str_split_fixed(somatic.complete.wf.zero.nic.outermost.lor.df$gamma.500.level, ' regions', n = 2)[, 1]

#### Construct iPSC data.frames ####
ipsc.complete.wf.all.lor.df <- rbind(ADS_iPSC.complete.wf.all.lor.df, IMR90_iPSC.complete.wf.all.lor.df, FF_iPSC_6.9.complete.wf.all.lor.df, FF_iPSC_19.7.complete.wf.all.lor.df, FF_iPSC_19.11.complete.wf.all.lor.df)
ipsc.complete.wf.all.lor.df$Sample <- c(rep('ADS-iPSC', nrow(ADS_iPSC.complete.wf.all.lor.df)), rep('IMR90-iPSC', nrow(IMR90_iPSC.complete.wf.all.lor.df)), rep('FF-iPSC_6.9', nrow(FF_iPSC_6.9.complete.wf.all.lor.df)), rep('FF-iPSC_19.7', nrow(FF_iPSC_19.7.complete.wf.all.lor.df)), rep('FF-iPSC_19.11', nrow(FF_iPSC_19.11.complete.wf.all.lor.df)))
# Remove " regions" from gamma.500.level labels to preserve space in plot facet labels
ipsc.complete.wf.all.lor.df$gamma.500.level <- str_split_fixed(ipsc.complete.wf.all.lor.df$gamma.500.level, ' regions', n = 2)[, 1]

ipsc.complete.wf.zero.nic.all.lor.df <- rbind(ADS_iPSC.complete.wf.zero.nic.all.lor.df, IMR90_iPSC.complete.wf.zero.nic.all.lor.df, FF_iPSC_6.9.complete.wf.zero.nic.all.lor.df, FF_iPSC_19.7.complete.wf.zero.nic.all.lor.df, FF_iPSC_19.11.complete.wf.zero.nic.all.lor.df)
ipsc.complete.wf.zero.nic.all.lor.df$Sample <- c(rep('ADS-iPSC', nrow(ADS_iPSC.complete.wf.zero.nic.all.lor.df)), rep('IMR90-iPSC', nrow(IMR90_iPSC.complete.wf.zero.nic.all.lor.df)), rep('FF-iPSC_6.9', nrow(FF_iPSC_6.9.complete.wf.zero.nic.all.lor.df)), rep('FF-iPSC_19.7', nrow(FF_iPSC_19.7.complete.wf.zero.nic.all.lor.df)), rep('FF-iPSC_19.11', nrow(FF_iPSC_19.11.complete.wf.zero.nic.all.lor.df)))
# Remove " regions" from gamma.500.level labels to preserve space in plot facet labels
ipsc.complete.wf.zero.nic.all.lor.df$gamma.500.level <- str_split_fixed(ipsc.complete.wf.zero.nic.all.lor.df$gamma.500.level, ' regions', n = 2)[, 1]

ipsc.complete.wf.outermost.lor.df <- rbind(ADS_iPSC.complete.wf.outermost.lor.df, IMR90_iPSC.complete.wf.outermost.lor.df, FF_iPSC_6.9.complete.wf.outermost.lor.df, FF_iPSC_19.7.complete.wf.outermost.lor.df, FF_iPSC_19.11.complete.wf.outermost.lor.df)
ipsc.complete.wf.outermost.lor.df$Sample <- c(rep('ADS-iPSC', nrow(ADS_iPSC.complete.wf.outermost.lor.df)), rep('IMR90-iPSC', nrow(IMR90_iPSC.complete.wf.outermost.lor.df)), rep('FF-iPSC_6.9', nrow(FF_iPSC_6.9.complete.wf.outermost.lor.df)), rep('FF-iPSC_19.7', nrow(FF_iPSC_19.7.complete.wf.outermost.lor.df)), rep('FF-iPSC_19.11', nrow(FF_iPSC_19.11.complete.wf.outermost.lor.df)))
# Remove " regions" from gamma.500.level labels to preserve space in plot facet labels
ipsc.complete.wf.outermost.lor.df$gamma.500.level <- str_split_fixed(ipsc.complete.wf.outermost.lor.df$gamma.500.level, ' regions', n = 2)[, 1]

ipsc.complete.wf.zero.nic.outermost.lor.df <- rbind(ADS_iPSC.complete.wf.zero.nic.outermost.lor.df, IMR90_iPSC.complete.wf.zero.nic.outermost.lor.df, FF_iPSC_6.9.complete.wf.zero.nic.outermost.lor.df, FF_iPSC_19.7.complete.wf.zero.nic.outermost.lor.df, FF_iPSC_19.11.complete.wf.zero.nic.outermost.lor.df)
ipsc.complete.wf.zero.nic.outermost.lor.df$Sample <- c(rep('ADS-iPSC', nrow(ADS_iPSC.complete.wf.zero.nic.outermost.lor.df)), rep('IMR90-iPSC', nrow(IMR90_iPSC.complete.wf.zero.nic.outermost.lor.df)), rep('FF-iPSC_6.9', nrow(FF_iPSC_6.9.complete.wf.zero.nic.outermost.lor.df)), rep('FF-iPSC_19.7', nrow(FF_iPSC_19.7.complete.wf.zero.nic.outermost.lor.df)), rep('FF-iPSC_19.11', nrow(FF_iPSC_19.11.complete.wf.zero.nic.outermost.lor.df)))
# Remove " regions" from gamma.500.level labels to preserve space in plot facet labels
ipsc.complete.wf.zero.nic.outermost.lor.df$gamma.500.level <- str_split_fixed(ipsc.complete.wf.zero.nic.outermost.lor.df$gamma.500.level, ' regions', n = 2)[, 1]

#### Construct ES data.frames ####
es.complete.wf.all.lor.df <- rbind(H1_r1.complete.wf.all.lor.df, H1_r2.complete.wf.all.lor.df, H9.complete.wf.all.lor.df, H9_Laurent.complete.wf.all.lor.df, HSF1.complete.wf.all.lor.df)
es.complete.wf.all.lor.df$Sample <- c(rep('H1_r1', nrow(H1_r1.complete.wf.all.lor.df)), rep('H1_r2', nrow(H1_r2.complete.wf.all.lor.df)), rep('H9', nrow(H9.complete.wf.all.lor.df)), rep('H9_Laurent', nrow(H9_Laurent.complete.wf.all.lor.df)), rep('HSF1', nrow(HSF1.complete.wf.all.lor.df)))
# Remove " regions" from gamma.500.level labels to preserve space in plot facet labels
es.complete.wf.all.lor.df$gamma.500.level <- str_split_fixed(es.complete.wf.all.lor.df$gamma.500.level, ' regions', n = 2)[, 1]

es.complete.wf.outermost.lor.df <- rbind(H1_r1.complete.wf.outermost.lor.df, H1_r2.complete.wf.outermost.lor.df, H9.complete.wf.outermost.lor.df, H9_Laurent.complete.wf.outermost.lor.df, HSF1.complete.wf.outermost.lor.df)
es.complete.wf.outermost.lor.df$Sample <- c(rep('H1_r1', nrow(H1_r1.complete.wf.outermost.lor.df)), rep('H1_r2', nrow(H1_r2.complete.wf.outermost.lor.df)), rep('H9', nrow(H9.complete.wf.outermost.lor.df)), rep('H9_Laurent', nrow(H9_Laurent.complete.wf.outermost.lor.df)), rep('HSF1', nrow(HSF1.complete.wf.outermost.lor.df)))
# Remove " regions" from gamma.500.level labels to preserve space in plot facet labels
es.complete.wf.outermost.lor.df$gamma.500.level <- str_split_fixed(es.complete.wf.outermost.lor.df$gamma.500.level, ' regions', n = 2)[, 1]

es.complete.wf.zero.nic.all.lor.df <- rbind(H1_r1.complete.wf.zero.nic.all.lor.df, H1_r2.complete.wf.zero.nic.all.lor.df, H9.complete.wf.zero.nic.all.lor.df, H9_Laurent.complete.wf.zero.nic.all.lor.df, HSF1.complete.wf.zero.nic.all.lor.df)
es.complete.wf.zero.nic.all.lor.df$Sample <- c(rep('H1_r1', nrow(H1_r1.complete.wf.zero.nic.all.lor.df)), rep('H1_r2', nrow(H1_r2.complete.wf.zero.nic.all.lor.df)), rep('H9', nrow(H9.complete.wf.zero.nic.all.lor.df)), rep('H9_Laurent', nrow(H9_Laurent.complete.wf.zero.nic.all.lor.df)), rep('HSF1', nrow(HSF1.complete.wf.zero.nic.all.lor.df)))
# Remove " regions" from gamma.500.level labels to preserve space in plot facet labels
es.complete.wf.zero.nic.all.lor.df$gamma.500.level <- str_split_fixed(es.complete.wf.zero.nic.all.lor.df$gamma.500.level, ' regions', n = 2)[, 1]

es.complete.wf.zero.nic.outermost.lor.df <- rbind(H1_r1.complete.wf.zero.nic.outermost.lor.df, H1_r2.complete.wf.zero.nic.outermost.lor.df, H9.complete.wf.zero.nic.outermost.lor.df, H9_Laurent.complete.wf.zero.nic.outermost.lor.df, HSF1.complete.wf.zero.nic.outermost.lor.df)
es.complete.wf.zero.nic.outermost.lor.df$Sample <- c(rep('H1_r1', nrow(H1_r1.complete.wf.zero.nic.outermost.lor.df)), rep('H1_r2', nrow(H1_r2.complete.wf.zero.nic.outermost.lor.df)), rep('H9', nrow(H9.complete.wf.zero.nic.outermost.lor.df)), rep('H9_Laurent', nrow(H9_Laurent.complete.wf.zero.nic.outermost.lor.df)), rep('HSF1', nrow(HSF1.complete.wf.zero.nic.outermost.lor.df)))
# Remove " regions" from gamma.500.level labels to preserve space in plot facet labels
es.complete.wf.zero.nic.outermost.lor.df$gamma.500.level <- str_split_fixed(es.complete.wf.zero.nic.outermost.lor.df$gamma.500.level, ' regions', n = 2)[, 1]

### Construct IVD data.frames ####
ivd.complete.wf.all.lor.df <- rbind(H1_BMP4.complete.wf.all.lor.df, FF_iPSC_19.11_BMP4.complete.wf.all.lor.df)
ivd.complete.wf.all.lor.df$Sample <- c(rep('H1+BMP4', nrow(H1_BMP4.complete.wf.all.lor.df)), rep('FF+iPSC_19.11_BMP4', nrow(FF_iPSC_19.11_BMP4.complete.wf.all.lor.df)))
# Remove " regions" from gamma.500.level labels to privderve space in plot facet labels
ivd.complete.wf.all.lor.df$gamma.500.level <- str_split_fixed(ivd.complete.wf.all.lor.df$gamma.500.level, ' regions', n = 2)[, 1]

ivd.complete.wf.outermost.lor.df <- rbind(H1_BMP4.complete.wf.outermost.lor.df, FF_iPSC_19.11_BMP4.complete.wf.outermost.lor.df)
ivd.complete.wf.outermost.lor.df$Sample <- c(rep('H1+BMP4', nrow(H1_BMP4.complete.wf.outermost.lor.df)), rep('FF+iPSC_19.11_BMP4', nrow(FF_iPSC_19.11_BMP4.complete.wf.outermost.lor.df)))
# Remove " regions" from gamma.500.level labels to privderve space in plot facet labels
ivd.complete.wf.outermost.lor.df$gamma.500.level <- str_split_fixed(ivd.complete.wf.outermost.lor.df$gamma.500.level, ' regions', n = 2)[, 1]

ivd.complete.wf.zero.nic.all.lor.df <- rbind(H1_BMP4.complete.wf.zero.nic.all.lor.df, FF_iPSC_19.11_BMP4.complete.wf.zero.nic.all.lor.df)
ivd.complete.wf.zero.nic.all.lor.df$Sample <- c(rep('H1+BMP4', nrow(H1_BMP4.complete.wf.zero.nic.all.lor.df)), rep('FF+iPSC_19.11_BMP4', nrow(FF_iPSC_19.11_BMP4.complete.wf.zero.nic.all.lor.df)))
# Remove " regions" from gamma.500.level labels to privderve space in plot facet labels
ivd.complete.wf.zero.nic.all.lor.df$gamma.500.level <- str_split_fixed(ivd.complete.wf.zero.nic.all.lor.df$gamma.500.level, ' regions', n = 2)[, 1]

ivd.complete.wf.zero.nic.outermost.lor.df <- rbind(H1_BMP4.complete.wf.zero.nic.outermost.lor.df, FF_iPSC_19.11_BMP4.complete.wf.zero.nic.outermost.lor.df)
ivd.complete.wf.zero.nic.outermost.lor.df$Sample <- c(rep('H1+BMP4', nrow(H1_BMP4.complete.wf.zero.nic.outermost.lor.df)), rep('FF+iPSC_19.11_BMP4', nrow(FF_iPSC_19.11_BMP4.complete.wf.zero.nic.outermost.lor.df)))
# Remove " regions" from gamma.500.level labels to privderve space in plot facet labels
ivd.complete.wf.zero.nic.outermost.lor.df$gamma.500.level <- str_split_fixed(ivd.complete.wf.zero.nic.outermost.lor.df$gamma.500.level, ' regions', n = 2)[, 1]

#### Change working directory to Lister_data.3x3.plots ####
setwd('~/Lister_data.3x3.plots')

#### Plots - Paired-end focus (full x-axis displayed), single-end focus (truncated x-axis). Both plots have truncated y-axesso that these are consistent across cell-types in order to make them comparable. This is done using coord_cartesian(limits = xxx) rather than scale_y_continuous(limits = xxx) to ensure the data aren't censored outside the limits (see https://groups.google.com/forum/?fromgroups=#!topic/ggplot2-dev/JGal_9nRzsg for details) ####
tmp <- ggplot(data = somatic.complete.wf.all.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'Somatic cell-types: Aggregated within-fragment comethylation\npair.choice = all') + scale_x_continuous('Distance between CpGs (bp)') + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.somatic.all_pairs.paired-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)
tmp <- ggplot(data = somatic.complete.wf.all.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'Somatic cell-types: Aggregated within-fragment comethylation\npair.choice = all') + scale_x_continuous('Distance between CpGs (bp)', breaks = c(0, 20, 40, 60, 80)) + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(xlim = c(0, 90), ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.somatic.all_pairs.single-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)

tmp <- ggplot(data = somatic.complete.wf.outermost.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'Somatic cell-types: Aggregated within-fragment comethylation\npair.choice = outermost') + scale_x_continuous('Distance between CpGs (bp)') + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.somatic.outermost_pairs.paired-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)
tmp <- ggplot(data = somatic.complete.wf.outermost.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'Somatic cell-types: Aggregated within-fragment comethylation\npair.choice = outermost') + scale_x_continuous('Distance between CpGs (bp)', breaks = c(0, 20, 40, 60, 80)) + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(xlim = c(0, 90), ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.somatic.outermost_pairs.single-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)

tmp <- ggplot(data = somatic.complete.wf.zero.nic.all.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'Somatic cell-types: Aggregated within-fragment comethylation\npair.choice = all, NIC = 0') + scale_x_continuous('Distance between CpGs (bp)') + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.somatic.all_pairs.zero_NIC.paired-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)
tmp <- ggplot(data = somatic.complete.wf.zero.nic.all.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'Somatic cell-types: Aggregated within-fragment comethylation\npair.choice = all, NIC = 0') + scale_x_continuous('Distance between CpGs (bp)', breaks = c(0, 20, 40, 60, 80)) + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(xlim = c(0, 90), ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.somatic.all_pairs.zero_NIC.single-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)

tmp <- ggplot(data = somatic.complete.wf.zero.nic.outermost.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'Somatic cell-types: Aggregated within-fragment comethylation\npair.choice = outermost, NIC = 0') + scale_x_continuous('Distance between CpGs (bp)') + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.somatic.outermost_pairs.zero_NIC.paired-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)
tmp <- ggplot(data = somatic.complete.wf.zero.nic.outermost.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'Somatic cell-types: Aggregated within-fragment comethylation\npair.choice = outermost, NIC = 0') + scale_x_continuous('Distance between CpGs (bp)', breaks = c(0, 20, 40, 60, 80)) + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(xlim = c(0, 90), ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.somatic.outermost_pairs.zero_NIC.single-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)

tmp <- ggplot(data = ipsc.complete.wf.all.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'iPSC cell-types: Aggregated within-fragment comethylation\npair.choice = all') + scale_x_continuous('Distance between CpGs (bp)') + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.ipsc.all_pairs.paired-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)
tmp <- ggplot(data = ipsc.complete.wf.all.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'iPSC cell-types: Aggregated within-fragment comethylation\npair.choice = all') + scale_x_continuous('Distance between CpGs (bp)', breaks = c(0, 20, 40, 60, 80)) + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(xlim = c(0, 90), ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.ipsc.all_pairs.single-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)

tmp <- ggplot(data = ipsc.complete.wf.outermost.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'iPSC cell-types: Aggregated within-fragment comethylation\npair.choice = outermost') + scale_x_continuous('Distance between CpGs (bp)') + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.ipsc.outermost_pairs.paired-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)
tmp <- ggplot(data = ipsc.complete.wf.outermost.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'iPSC cell-types: Aggregated within-fragment comethylation\npair.choice = outermost') + scale_x_continuous('Distance between CpGs (bp)', breaks = c(0, 20, 40, 60, 80)) + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(xlim = c(0, 90), ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.ipsc.outermost_pairs.single-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)

tmp <- ggplot(data = ipsc.complete.wf.zero.nic.all.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'iPSC cell-types: Aggregated within-fragment comethylation\npair.choice = all, NIC = 0') + scale_x_continuous('Distance between CpGs (bp)') + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.ipsc.all_pairs.zero.nic.paired-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)
tmp <- ggplot(data = ipsc.complete.wf.zero.nic.all.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'iPSC cell-types: Aggregated within-fragment comethylation\npair.choice = all, NIC = 0') + scale_x_continuous('Distance between CpGs (bp)', breaks = c(0, 20, 40, 60, 80)) + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(xlim = c(0, 90), ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.ipsc.all_pairs.zero.nic.single-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)

tmp <- ggplot(data = ipsc.complete.wf.zero.nic.outermost.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'iPSC cell-types: Aggregated within-fragment comethylation\npair.choice = outermost, NIC = 0') + scale_x_continuous('Distance between CpGs (bp)') + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.ipsc.outermost_pairs.zero.nic.paired-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)
tmp <- ggplot(data = ipsc.complete.wf.zero.nic.outermost.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'iPSC cell-types: Aggregated within-fragment comethylation\npair.choice = outermost, NIC = 0') + scale_x_continuous('Distance between CpGs (bp)', breaks = c(0, 20, 40, 60, 80)) + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(xlim = c(0, 90), ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.ipsc.outermost_pairs.zero.nic.single-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)

tmp <- ggplot(data = es.complete.wf.all.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'ES cell-types: Aggregated within-fragment comethylation\npair.choice = all') + scale_x_continuous('Distance between CpGs (bp)') + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.es.all_pairs.paired-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)
tmp <- ggplot(data = es.complete.wf.all.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'ES cell-types: Aggregated within-fragment comethylation\npair.choice = all') + scale_x_continuous('Distance between CpGs (bp)', breaks = c(0, 20, 40, 60, 80)) + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(xlim = c(0, 90), ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.es.all_pairs.single-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)

tmp <- ggplot(data = es.complete.wf.outermost.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'ES cell-types: Aggregated within-fragment comethylation\npair.choice = outermost') + scale_x_continuous('Distance between CpGs (bp)') + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.es.outermost_pairs.paired-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)
tmp <- ggplot(data = es.complete.wf.outermost.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'ES cell-types: Aggregated within-fragment comethylation\npair.choice = outermost') + scale_x_continuous('Distance between CpGs (bp)', breaks = c(0, 20, 40, 60, 80)) + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(xlim = c(0, 90), ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.es.outermost_pairs.single-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)

tmp <- ggplot(data = es.complete.wf.zero.nic.all.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'ES cell-types: Aggregated within-fragment comethylation\npair.choice = all, NIC = 0') + scale_x_continuous('Distance between CpGs (bp)') + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.es.all_pairs.zero.nic.paired-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)
tmp <- ggplot(data = es.complete.wf.zero.nic.all.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'ES cell-types: Aggregated within-fragment comethylation\npair.choice = all, NIC = 0') + scale_x_continuous('Distance between CpGs (bp)', breaks = c(0, 20, 40, 60, 80)) + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(xlim = c(0, 90), ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.es.all_pairs.zero.nic.single-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)

tmp <- ggplot(data = es.complete.wf.zero.nic.outermost.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'ES cell-types: Aggregated within-fragment comethylation\npair.choice = outermost, NIC = 0') + scale_x_continuous('Distance between CpGs (bp)') + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.es.outermost_pairs.zero.nic.paired-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)
tmp <- ggplot(data = es.complete.wf.zero.nic.outermost.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'ES cell-types: Aggregated within-fragment comethylation\npair.choice = outermost, NIC = 0') + scale_x_continuous('Distance between CpGs (bp)', breaks = c(0, 20, 40, 60, 80)) + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(xlim = c(0, 90), ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.es.outermost_pairs.zero.nic.single-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)

tmp <- ggplot(data = ivd.complete.wf.all.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'IVD cell-types: Aggregated within-fragment comethylation\npair.choice = all') + scale_x_continuous('Distance between CpGs (bp)') + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.ivd.all_pairs.paired-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)
tmp <- ggplot(data = ivd.complete.wf.all.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'IVD cell-types: Aggregated within-fragment comethylation\npair.choice = all') + scale_x_continuous('Distance between CpGs (bp)', breaks = c(0, 20, 40, 60, 80)) + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(xlim = c(0, 90), ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.ivd.all_pairs.single-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)

tmp <- ggplot(data = ivd.complete.wf.outermost.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'IVD cell-types: Aggregated within-fragment comethylation\npair.choice = outermost') + scale_x_continuous('Distance between CpGs (bp)') + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.ivd.outermost_pairs.paired-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)
tmp <- ggplot(data = ivd.complete.wf.outermost.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'IVD cell-types: Aggregated within-fragment comethylation\npair.choice = outermost') + scale_x_continuous('Distance between CpGs (bp)', breaks = c(0, 20, 40, 60, 80)) + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(xlim = c(0, 90), ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.ivd.outermost_pairs.single-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)

tmp <- ggplot(data = ivd.complete.wf.zero.nic.all.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'IVD cell-types: Aggregated within-fragment comethylation\npair.choice = all, NIC = 0') + scale_x_continuous('Distance between CpGs (bp)') + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.ivd.all_pairs.zero_NIC.paired-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)
tmp <- ggplot(data = ivd.complete.wf.zero.nic.all.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'IVD cell-types: Aggregated within-fragment comethylation\npair.choice = all, NIC = 0') + scale_x_continuous('Distance between CpGs (bp)', breaks = c(0, 20, 40, 60, 80)) + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(xlim = c(0, 90), ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.ivd.all_pairs.zero_NIC.single-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)

tmp <- ggplot(data = ivd.complete.wf.zero.nic.outermost.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'IVD cell-types: Aggregated within-fragment comethylation\npair.choice = outermost, NIC = 0') + scale_x_continuous('Distance between CpGs (bp)') + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.ivd.outermost_pairs.zero_NIC.paired-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)
tmp <- ggplot(data = ivd.complete.wf.zero.nic.outermost.lor.df, aes(x = IPD, y = lor, colour = Sample)) + facet_grid(CGI ~ gamma.500.level) + geom_line(size = 2) + presentation.theme + ggtitle(label = 'IVD cell-types: Aggregated within-fragment comethylation\npair.choice = outermost, NIC = 0') + scale_x_continuous('Distance between CpGs (bp)', breaks = c(0, 20, 40, 60, 80)) + scale_y_continuous('Log odds ratio', breaks = 0:4) + coord_cartesian(xlim = c(0, 90), ylim = c(-1, 4.5))
ggsave(filename = 'Lister_data.ivd.outermost_pairs.zero_NIC.single-end_focus.3x3.pdf', plot = tmp, width = 16.9, height = 10.5)
rm(tmp)

#### Finished ####


