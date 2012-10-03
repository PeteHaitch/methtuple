# Peter Hickey
# 02/10/2012
# Aggregated-by-distance within-fragment comethylation plotted by cell-type

#### TODOs ####
# Need to stratify by tertiles of gamma before doing these plots

#### Load libraries ####
library(ggplot2)

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
rm(ADS_adipose.env)
gc()

ADS_iPSC.env <- new.env()
load('Lister_2011_BS-seq_data/ADS-iPSC/Downloaded_mapped_reads/WF/ADS-iPSC_comethylation.RData', envir = ADS_iPSC.env)
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
rm(ADS_iPSC.env)
gc()

FF.env <- new.env()
load('Lister_2011_BS-seq_data/FF/Downloaded_mapped_reads/WF/FF_comethylation.RData', envir =  FF.env)
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
rm(FF.env)
gc()

FF_iPSC_19.11.env <- new.env()
load('Lister_2011_BS-seq_data/FF-iPSC_19.11/Downloaded_mapped_reads/WF/FF-iPSC_19.11_comethylation.RData', envir =  FF_iPSC_19.11.env)
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
rm(FF_iPSC_19.11.env)
gc()

FF_iPSC_19.11_BMP4.env <- new.env()
load('Lister_2011_BS-seq_data/FF-iPSC_19.11+BMP4/Downloaded_mapped_reads/WF/FF-iPSC_19.11+BMP4_comethylation.RData', envir =  FF_iPSC_19.11_BMP4.env)
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
rm(FF_iPSC_19.11_BMP4.env)
gc()

FF_iPSC_19.7.env <- new.env()
load('Lister_2011_BS-seq_data/FF-iPSC_19.7/Downloaded_mapped_reads/WF/FF-iPSC_19.7_comethylation.RData', envir =  FF_iPSC_19.7.env)
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
rm(IMR90_iPSC.env)
gc()

save.image('bloop.RData')

somatic.all.lor <- rbind(IMR90_r1.all.lor, IMR90_r2.all.lor, FF.all.lor, ADS.all.lor, ADS_adipose.all.lor)
somatic.all.lor$Sample <- c(rep('IMR90_r1', nrow(IMR90_r1.all.lor)), rep('IMR90_r2', nrow(IMR90_r2.all.lor)), rep('FF', nrow(FF.all.lor)), rep('ADS', nrow(ADS.all.lor)), rep('ADS-adipose', nrow(ADS_adipose.all.lor)))
somatic$CGI <- ifelse(somatic$CGI == 'CGI', 'CpG island', 'Non CpG island')