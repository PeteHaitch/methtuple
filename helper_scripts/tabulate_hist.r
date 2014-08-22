#### DESCRIPTION ####
# Peter Hickey (peter.hickey@gmai.com)
# 05/02/2014
# Merge .hist files created by run_methtuple.sh helper script

#### Parse command line argument to set the working directory ####
wd <- commandArgs(TRUE)[1]
setwd(wd)

#### Get filenames of .hist files and read this files into a list ####
hist_filenames <- list.files(pattern = "\\.hist$")
hist_files <- lapply(X = hist_filenames, FUN = function(x){read.table(x, header = T, sep = '\t', as.is = TRUE)})

#### Identify unique (sorted) 'n'; sum all 'count' values with n = 'n'; construct 'd', a data.frame of resulting aggregate 'n' and 'count' ####
n <- unique(sort(do.call(c, lapply(hist_files, FUN = function(hist_file){hist_file[, 'n']}))))
count <- sapply(n, FUN = function(nn, hist_files){
  sum(as.numeric(unlist(lapply(hist_files, FUN = function(hist_file, nn){
    hist_file[hist_file[, 'n'] == nn, 'count']
    }, nn = nn))))
  }, hist_files = hist_files)
d <- data.frame(n = n, count = count)

#### Write 'd' to file ####
# Re-construct the sample name, which basically means removing _<CHROM> part from from the first part of the hist_filenames
sample_name <- paste0(head(strsplit(strsplit(hist_filenames[[1]], '\\.')[[1]][1], '_')[[1]], -1), collapse = '_')
# The middle part of the filename is "<methylationType>_per_read
middle_part <- strsplit(hist_filenames[[1]], '\\.')[[1]][[2]]
# Construct the complete out filename
out_filename <- paste(sample_name, middle_part, 'hist', sep = '.')
write.table(d, file = out_filename, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
