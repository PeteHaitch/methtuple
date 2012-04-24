#!/usr/bin/python
import re
import random
from numpy import array
import argparse
import sys
import csv
import pysam

# Input SAM/BAM must be query-name-sorted
# Input SAM/BAM must have been corrected with correct_Bismark_PE_SAM.py. correct_Bismark_PE_SAM.py only works for directional libraries.
# A read-pair is comprised of read1 and read2
# If a read1 and read2 overlap (i.e. the DNA fragment was shorter then the sum of the read lenghts) I exclude the overlapping positions from the 3' end of read2 since read2 is generally of lower quality. A better approach would be to trim bases from either read1 or read2 depending on which read had the lower quality 'overlap-bases'. The Epigenomic Roadmap recommend that "Paired-reads that have partial overlap in genome coverage should be trimmed from the 3' so as to avoid treating sequence derived from multiple passes of the same genomic DNA fragment as independent data points." (http://www.roadmapepigenomics.org/files/protocols/data/dna-methylation/MethylC-SeqStandards_FINAL.pdf). A further improvement would be a check that the overlap reports a consistent methylation call. A further complication is addressed by bad_overlap_pairs.
# TLEN is defined as the rightmost position of read2 (read1) minus the leftmost position of read1 (read2) plus 1 for OT (OB) aligned reads. Thus if TLEN < read1.alen + read2.alen the two reads overlap.

## TODO: Does read1 always appear before read2 in a query-name-sorted BAM?
## TODO: Output all CpG-pairs in a read. Pairs are created using the closest, rightmost CpG; i.e. a read with 2 CpGs has 1 pair, 3 CpGs => 2 pairs, 4 CpGs => 3 pairs, ..., n CpGs => (n-1) pairs
## TOOD: Check if MarkDuplicates marks both reads of a read-pair as duplicates
## TODO: Add a count of the number of CpGs between the CpG-pair to the output for each CpG pair
## TODO: Add flags to skip 'non-primary alignment' reads, 'failed reads', etc.
## TODO: Bismark's TLEN computation is incorrect for overlapping reads from a read-pair. Thus it should not be used for inferring that reads overlap. Rewrite any features that rely on inferring read overlaps so as not to rely on the TLEN value. A bug report has been submitted to the developer's of Bismark.
## TODO: Fix bug in calculation of overlapping reads - reads overlap if abs(TLEN) < read1.alen + read2.alen
## TODO: Fix bug when ignoring reads that have "bad" overlap.
## TODO: Add parameter to specify minimum and maximum insert sizes for a read-pair to be considered
## TODO: Write bad_overlap_reads to a separate file.
## TODO: Make sure that Case1C and Case2B point to the C in the CpG and not the G
## TODO: Test-case with reads that contain an overlap.
## TODO: Add counters for each of the cases, CpGs_per_read, CpG_positions, methylated_CpG_positions, unmethylated_CpG_positions, bad_overlap_pairs
## TODO: Add --ignore5 and --ignore3 options
## TODO: Add check if --ignoreDuplicates is set and ignore duplicates if it is set

### Command line passer
parser = argparse.ArgumentParser(description='Extract the methylation calls for a CpG pair from reads that overlap multiple CpGs from a Bismark SAM/BAM file. The output file contains the positions of the CpG pair using 1-based co-ordinates on the forward strand of the cytosine in each CpG. If a read overlaps more than two CpGs there are several ways to construct the pairs (see the --pairChoice argument). The output of this file can be used for analysing comethylation along a read.')
parser.add_argument('SAM', metavar = 'SAM',
                  help='The path to the Bismark SAM/BAM file')
parser.add_argument('output', type=argparse.FileType('w'), nargs=1,
                   help='The output filename')
args = parser.parse_args()

### Create (possibly redundant) pointers to files from command line arguments
OUT = args.output[0]

### Open SAM/BAM file
f = pysam.Samfile(args.SAM) 

### Variable initialisations
CpG_pattern = re.compile(r"[Zz]")
methylated_CpG_pattern = re.compile(r"[Z]")
unmethylated_CpG_pattern = re.compile(r"[z]")
ignore2 = 0 # The number of bases of read2 to ignore due to overlap between read1 and read2
bad_overlap_pairs = 0 # The number of read-pairs where the leftmost position of read2 (read1) is to the left of the leftmost position of read1 (read2) for reads informative for the OT (OB) strand. Such reads are ignored when computing single-strand comethylation as they can violate the assumption that posL < posR.

### Function definitions
## getMate() returns a pointer to read2 of a properly-paired read-pair
def getMate(read1Pointer): # read1Pointer is a pointer in the BAM file to the first read in a properly-paired read-pair
    return read2Pointer

## SE_SAM2MS extracts and summarises the methylation string (MS) information for a read mapped as single-end data (i.e. with the appropriate SAM flag set)
def SE_SAM2MS(read):
    print XM
    return pointerToNextReadPair

## PE_SAM2MS extracts and summarises the methylation string (MS) information for a read-pair mapped as paired-end data (i.e. with the appropriate SAM flag set).
def PE_SAM2MS(read, mate):
    mate = getMate(read1Pointer)
    print XM
    return pointerToNextReadPair

## tabWriter writes a tab-separated output file
tabWriter = csv.writer(OUT, delimiter='\t', quotechar=' ', quoting=csv.QUOTE_MINIMAL)

### Loop over the BAM file line-by-line (i.e. alignedRead-by-alignedRead)
for linePointer in f: # Pseudo-code; will not work as intended.
    read = getRead(linePointer)
    if read.is_pair:
        PE_SAM2MS(read)
    else:
        SE_SAM2MS(read)
        
        



f.close()
OUT.close()
