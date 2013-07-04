#### Class definition ####
# Peter Hickey
# 30/01/2013
# WFComethylation class for storing *.wf files
# What I want: (1) Based on GRanges; (2) sampleName included; (3) nTuple included; (4) slots for n-tuple-specific and aggregated log odds ratios

#### contingencyTable ####
# TODO: Either define a method 'table' for WFComethylation objects or create a contingencyTable class which works on WFComethylation objects 

#### Definition ####
# aggregate log odds ratios (LOR) for a single sample and strata
# e.g. for sample 'Cancer' within CpG islands 
setClass("aggregate_LOR",
         representation(
           sampleName = "character",
           strata = 'character'
           contingencyTable = "array"
           methylationType = "character"
         ))

#### Validity checker ####
setValidity("aggregate_LOR", function(object) {
  msg <- NULL
  ### UP TO HERE ###
  if (is.null(names(object@contingencyTable)))
    msg <- c(msg, "'contingenctTable' must have names")
  nNTuples <- length(object@gr)
  if(nNTuples != nrow(object@ntuples) || nNTuples != nrow(object@positions))
    msg <- c(msg, "length of the 'gr' slot is not equal to the nrows of the 'ntuples' or 'positions' slots")
  if (any(object@ntuples < 0))
    msg <- c(msg, "the 'ntuples' slot has negative entries")
  if (!is.null(object@NIML))
    if ((length(object@NIML)) != nrow(object@positions))
      msg <- c(msg, "length of the 'NIML' slot is not equal to nrow of the 'positions' slot")
  if (length(object@sampleName) != 1)
    msg <- c(msg, "length of 'sampleName' is not equal to 1")
  if (!(object@methylationType %in% c("CpG", "CHG", "CHH")))
    msg <- c(msg, "'methylationType' must be one of 'CpG', 'CHG' or 'CHH'")
  if(is.null(msg)) TRUE else msg
})