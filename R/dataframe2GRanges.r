# Written by Kasper Hansen
# https://stat.ethz.ch/pipermail/bioconductor/2010-September/035341.html

data.frame2GRanges <- function(object, keepColumns = TRUE) {
    stopifnot(class(object) == "data.frame")
    stopifnot(all(c("chr", "start", "end") %in% names(object)))
    if("strand" %in% names(object)) {
        if(is.numeric(object$strand)) {
            strand <- ifelse(object$strand == 1, "+", "*")
            strand[object$strand == -1] <- "-"
            object$strand <- strand
        }
        gr <- GRanges(seqnames = object$chr,
                      ranges = IRanges(start = object$start, end = object$end),
                      strand = object$strand)
    } else {
        gr <- GRanges(seqnames = object$chr,
                      ranges = IRanges(start = object$start, end = object$end))
    }
    if(keepColumns) {
        dt <- as(object[, setdiff(names(object), c("chr", "start",
"end", "strand"))],
                     "DataFrame")
        elementMetadata(gr) <- dt
    }
    names(gr) <- rownames(object)
    gr
}


### Something interesting Kasper wrote in a lecture (http://www.bioconductor.org/help/course-materials/2011/CSAMA/Tuesday/Morning%20Talks/IRangesLecture.pdf)
# "Usecase II
# Suppose we have a set of DMRs (think genomic regions) and a set of
# CpG Islands and we want to nd all DMRs within 10kb of a CpG Island.
# dmrs: GRanges
# islands: GRanges
# > findOverlaps(dmrs, resize(islands,
# + width = 20000 + width(islands), fix = "center"))
# (watch out for strand)"