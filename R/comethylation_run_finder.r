# Peter Hickey
# 09/10/2012
# Find runs of comethylation to try to discover interesting regions of the genome

#' Find the longest run of comethylation > threshold from a wf.zero.nic.all.gr object
#' 
#' @param x is an wf.zero.nic.all.gr object
#' @param threshold is the value of the log odds ratio we require to be greater then (less than) to define a run
#' @param positive is TRUE if we want to be greater than the threshold, FALSE if we want to be less than the threshold
#' 
#' @note Based on http://stackoverflow.com/questions/4655848/calculating-a-consecutive-streak-in-data
findMaxRunLengthOfComethylation <- function(x, threshold, greater.than = TRUE, na.rm = TRUE) {
  if(greater.than) {
    rles <- rle(values(x)$lor > threshold)
    max(rles$lengths[rles$values], na.rm = na.rm)
  } else{
    rles <- rle(values(x)$lor < threshold)
    max(rles$lengths[!rles$values], na.rm = na.rm)
  }
}