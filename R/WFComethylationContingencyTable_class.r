#### Class definition ####
# Peter Hickey
# 1/02/2013
# WFComethylationContingencyTables class for storing contingency tables based on a WFComethylation object

#### contingencyTable ####
setClass("WFComethylationContingencyTables",
         representation(
           sampleName = "character",
           stratificationName = "character",
           strata = "list",           
           methylationType = "character"
         ))