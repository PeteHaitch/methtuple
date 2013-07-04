#### Utility definitions for working with WFComethylation class objects ####
# Peter Hickey
# 30/01/2013

#' @param WFComethylation is a WFComethylation class object of the sample
#' @param methylation_loci_bed_file is the path to a bed file of methylation loci in the sample or its reference genome.
#' 
#' @return An updated version of the WFComethylation object with a vector of the number of intervening methylation loci (NIML) added to the 'NIML' slot. WARNING: If the NIML slot already exists it will be overwritten in the returned object.
#' @export
add_NIML <- function(WFComethylation, methylation_loci_bed_file){
  stopifnot(class(WFComethylation) == "WFComethylation")
  ## Import the bed file of CpGs
  cat('Reading in', methylation_loci_bed_file, '...\n')
  reading_time <- system.time(methylation_loci <- import(con = methylation_loci_bed_file, asRangedData = FALSE))[3]    
  cat('Finished reading ', methylation_loci_bed_file, 'in', reading_time, 'seconds.\n')
  cat('Counting NIML for each n-tuple...\n')
  ## Count the NIML for each methylation loci n-tuple
  n <- ncol(WFComethylation) # The 'n' in n-tuple
  counting_time <- system.time(WFComethylation@NIML <- countOverlaps(WFComethylation, methylation_loci) - n)[3] # All n-tuples overlap n methylation loci by definition, therefore need to subtract n from each total.
  cat('Finished counting NIML for each n-tuple in', counting_time, 'seconds.\n')
  return(WFComethylation)
}

# TODO: Function to compute n-tuple-specific log odds ratios
# TODO: Function to compute aggregated log odds ratios
