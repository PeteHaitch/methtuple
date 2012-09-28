# How to load RData to a specific environment, e.g.
H1_r1.env <- new.env()
load('Lister_2009_BS-seq_data/H1/r1/Downloaded_mapped_reads/WF/H1_r1_comethylation.RData', envir=H1_r2.env)
H1_r2.env <- new.env()
load('Lister_2009_BS-seq_data/H1/r2/Downloaded_mapped_reads/WF/H1_r2_comethylation.RData', envir=H1_r2.env)

# And to access the loaded objects
H1_r1.env$all.lor
H1_r2.env$all.lor

