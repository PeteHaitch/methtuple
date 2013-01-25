##### DESCRIPTION #####
# Peter Hickey 
# 24/01/2013
# The following is all taken from bsseq v0.7.3
setClassUnion("matrixOrNULL", c("matrix", "NULL"))

setClass("hasGRanges",
         representation(gr = "GRanges"))

setMethod("seqnames", signature(x = "hasGRanges"), function(x) {
    seqnames(x@gr)
})
setReplaceMethod("seqnames", "hasGRanges", function(x, value) {
    gr <- granges(x)
    seqnames(gr) <- value
    x@gr <- gr
    x
})

setMethod("seqlevels", signature(x = "hasGRanges"), function(x) {
    seqlevels(x@gr)
})
setReplaceMethod("seqlevels", "hasGRanges", function(x, value) {
    gr <- granges(x)
    seqlevels(gr) <- value
    x@gr <- gr
    x
})

setMethod("seqlengths", signature(x = "hasGRanges"), function(x) {
    seqlengths(x@gr)
})
setReplaceMethod("seqlengths", "hasGRanges", function(x, value) {
    gr <- granges(x)
    seqlengths(gr) <- value
    x@gr <- gr
    x
})

setMethod("granges", signature(x = "hasGRanges"),
          function(x) x@gr)

## FIXME: might want a granges replacement function

setMethod("start", "hasGRanges", function(x, ...) {
    start(x@gr, ...)
})
setReplaceMethod("start", "hasGRanges", function(x, check = TRUE, value) {
    gr <- granges(x)
    start(gr, check = check) <- value
    x@gr <- gr
    x
})

setMethod("end", "hasGRanges", function(x, ...) {
    end(x@gr, ...)
})
setReplaceMethod("end", "hasGRanges", function(x, check = TRUE, value) {
    gr <- granges(x)
    end(gr, check = check) <- value
    x@gr <- gr
    x
})

setMethod("width", "hasGRanges", function(x) {
    width(x@gr)
})
setReplaceMethod("width", "hasGRanges", function(x, check = TRUE, value) {
    gr <- granges(x)
    width(gr, check = check) <- value
    x@gr <- gr
    x
})

setMethod("strand", "hasGRanges", function(x) {
    strand(x@gr)
})
setReplaceMethod("strand", "hasGRanges", function(x, value) {
    gr <- granges(x)
    strand(gr) <- value
    x@gr <- gr
    x
})

setMethod("length", "hasGRanges", function(x) length(x@gr))



setMethod("[", "hasGRanges", function(x, i, ...) {
    if(missing(i))
        stop("need [i] for subsetting")
    if(missing(i))
        return(x)
    x@gr <- x@gr[i]
    x
})

#### End of things taken from bsseq ####

