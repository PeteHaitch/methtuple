#!/usr/bin/python
import re
import random
from numpy import array
import argparse
import sys
import csv
import pysam

### Notes and comments
##############################################################################################################################################################################################
# Input SAM/BAM must be query-name-sorted.
# Input SAM/BAM must have been corrected with correct_Bismark_PE_SAM.py. correct_Bismark_PE_SAM.py only works for directional libraries.
# A read-pair is comprised of read1 and read2
# Read1 appears before read2 in a query-name-sorted BAM
# Picard's MarkDuplicates flags paired-end reads as duplicates if they have identical 5' position of the leftmost read and identical TLENs.
# If a read1 and read2 overlap (i.e. the DNA fragment was shorter then the sum of the read lenghts) I exclude the overlapping positions from the 3' end of read2 since read2 is generally of lower quality. A better approach would be to trim bases from either read1 or read2 depending on which read had the lower quality 'overlap-bases'. The Epigenomic Roadmap recommend that "Paired-reads that have partial overlap in genome coverage should be trimmed from the 3' so as to avoid treating sequence derived from multiple passes of the same genomic DNA fragment as independent data points." (http://www.roadmapepigenomics.org/files/protocols/data/dna-methylation/MethylC-SeqStandards_FINAL.pdf). A further improvement would be a check that the overlap reports a consistent methylation call. A further complication is addressed by bad_overlap_pairs.
# TLEN is defined as the rightmost position of read2 (read1) minus the leftmost position of read1 (read2) plus 1 for OT (OB) aligned reads. Thus if TLEN < read1.alen + read2.alen the two reads overlap.
##############################################################################################################################################################################################

##############################################################################################################################################################################################
# WARNING: SCRIPT IS WRONG BECAUSE STRAND INFORMATION IS INCORRECTLY EXTRACTED
##############################################################################################################################################################################################

### TODOs ###
##############################################################################################################################################################################################
## TODO: Extract strand information from XS tag, rather than relying on read orientation (which is the incorrect way to infer strand information).
## TODO: Output all CpG-pairs in a read. Pairs are created using the closest, rightmost CpG; i.e. a read with 2 CpGs has 1 pair, 3 CpGs => 2 pairs, 4 CpGs => 3 pairs, ..., n CpGs => (n-1) pairs
## TODO: Check if MarkDuplicates marks both reads of a read-pair as duplicates
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
##############################################################################################################################################################################################

### Command line passer ###
##############################################################################################################################################################################################
parser = argparse.ArgumentParser(description='Extract the methylation calls for a CpG pair from reads that overlap multiple CpGs from a Bismark SAM/BAM file. The output file contains the positions of the CpG pair using 1-based co-ordinates on the forward strand of the cytosine in each CpG. If a read overlaps more than two CpGs there are several ways to construct the pairs (see the --pairChoice argument). The output of this file can be used for analysing comethylation along a read.')
parser.add_argument('BAM', metavar = 'BAM',
                  help='The path to the Bismark SAM/BAM file')
parser.add_argument('output', type=argparse.FileType('w'), nargs=1,
                   help='The output filename')
parser.add_argument('--ignoreDuplicates', action='store_true',help='Ignore reads that have been flagged as PCR duplicates by Picard\'s MarkDuplicates function')
args = parser.parse_args()
##############################################################################################################################################################################################

### Open SAM/BAM file and create (possibly redundant) pointers to output file from command line arguments ###
##############################################################################################################################################################################################
BAM = pysam.Samfile(args.BAM)
OUT = args.output[0]
##############################################################################################################################################################################################

### Variable initialisations ###
##############################################################################################################################################################################################
CpG_pattern = re.compile(r"[Zz]")
methylated_CpG_pattern = re.compile(r"[Z]")
unmethylated_CpG_pattern = re.compile(r"[z]")
bad_overlap_pairs = 0 # The number of read-pairs with a bad overlap; see PE_SAM2MS() for the definition of a bad overlap.
##############################################################################################################################################################################################

### Function definitions ###
##############################################################################################################################################################################################
## SE_SAM2MS extracts and summarises the methylation string (MS) information for a read mapped as single-end data (i.e. with the appropriate SAM flag set)
def SE_SAM2MS(read):
    return XM

## PE_SAM2MS extracts and summarises the methylation string (MS) information for a read-pair mapped as paired-end data (i.e. with the appropriate SAM flag set).
# Returns one of three integer values: -1, 0, 1. The XM-tag information is written to the output file (OUT) as a side-effect of calling this function.
# A return value of 1 means the XM information was successfully summarised and written to the output file (OUT). This includes the case where the read-pair did not overlap > 1 CpG and thus the read-pair is uninformative for CpG comethylation; nothing is written to the output file (OUT) in this latter case.
# A return value of 0 means there was a "bad overlap" between read1 and read2. A bad overlap for a read-pair that is informative for the OT strand is one where the leftmost (3') position of read2 is to the left of the leftmost (5') position of read1. Conversely, a bad overlap for a read-pair that is informative for the OB strand is one where the leftmost (3') position of read1 is to the left of the leftmost (5') position of read2. Such read-pairs are excluded from the analysis since it is more difficult to extract the XM information from these read-pairs, and hence no XM information was written to the output file (OUT).
# A return value of -1 means there was an error in extracting the XM information. This could be due to a an unexpected orientation of the read-pair or an unexpected distribution of CpGs in the read-pair.

def PE_SAM2MS(read1, read2):
    chrom = BAM.getrname(read1.tid)
    ignore2 = 0 # Initialise ignore2 to zero. ignore2 = the number of bases of read2 to ignore due to overlap between read1 and read2 and is computed later in this function.
    
    # Store the XM tag for each read. Tags are stored as a Python list of 2-tuples [("TAG_NAME", "TAG_VALUE"), ...]. This assumes there is one, and one only, XM tag per read
    read1XM = [tag[1] for tag in read1.tags if tag[0]=='XM'][0] 
    read2XM = [tag[1] for tag in read2.tags if tag[0]=='XM'][0]
    
    # Case1: read-pair informative for OT strand, read-pair orientation is +/- (read1/read2)
    if (not read1.is_reverse and read2.is_reverse):
        
        # If read1.tlen < 0, then the read-pair has a bad overlap
        if read1.tlen < 0:
            print 'WARNING: Ignoring read-pair ', read1.qname, ' due to bad overlap (=', read1.tlen, ')'
            return 0
        
        # start1 <= start2 by definition for a properly paired OT read
        start1 = read1.pos + 1 # read1.pos is 0-based in accordance with BAM file specificiations, regardless of whether the input file is SAM or BAM. I convert this to a 1-based position.
        start2 = read2.pos + 1 # read2.pos is 0-based in accordance with BAM file specificiations, regardless of whether the input file is SAM or BAM. I convert this to a 1-based position.
        strand = "+"

        # If the two reads overlap, i.e the rightmost position of read1 is greater than or equal to the leftmost position of read2, drop the offending positions from the leftmost (3') position of read2
        if read1.positions[-1] >= read2.positions[0]: 
            ignore2 = int(read1.positions[-1] - read2.positions[0] + 1) # In this case, ignore2 is the number of bases to ignore from the 3' end of read2.
            read2XM = read2XM[ignore2:] # Drop the ignored positions from read2's XM tag
            
        # Identify CpGs in read1 and read2
        CpG_index1 = [m.start() for m in re.finditer(CpG_pattern, read1XM)]
        CpG_index2 = [m.start() for m in re.finditer(CpG_pattern, read2XM)]
        
        # CaseA: > 0 CpGs in both reads of pair    
        if (len(CpG_index1) > 0 and len(CpG_index2) > 0):
            CpGL = CpG_index1[0] # Leftmost CpG in read1
            CpGR = CpG_index2[-1] # Rightmost CpG in read2
            positionL = start1 + CpGL
            positionR = start2 + CpGR + ignore2 # +ignore2 to deal with ignored bases
            output=[chrom, positionL, positionR, read1XM[CpGL], read2XM[CpGR], strand] # Output is left-to-right, thus read1XM appears before read2XM for OT read-pairs.
            if positionL > positionR:
                print "ERROR: Case1A posL > posR for read-pair ", read1.qname," with output ", output
                return -1
            else:
                tabWriter.writerow(output)
                return 1
        # CaseB: > 1 CpG in read1 but 0 CpGs in read2
        elif (len(CpG_index1) > 1 and len(CpG_index2) == 0):
            CpGL = CpG_index1[0]
            CpGR = CpG_index1[-1]
            positionL = start1 + CpGL
            positionR = start1 + CpGR
            output=[chrom, positionL, positionR, read1XM[CpGL], read1XM[CpGR], strand]
            if positionL > positionR:
                print "ERROR: Case1B posL > posR for read ", read1.qname," with output ", output
                return -1
            else:
                tabWriter.writerow(output)
                return 1
        # CaseC: 0 CpGs in read1 but > 1 CpG in read2
        elif (len(CpG_index1) == 0 and len(CpG_index2) > 1):
            CpGL = CpG_index2[0]
            CpGR = CpG_index2[-1]
            positionL = start2 + CpGL
            positionR = start2 + CpGR
            output=[chrom, positionL, positionR, read2XM[CpGL], read2XM[CpGR], strand]
            if positionL > positionR:
                print "ERROR: Case1C posL > posR for read ", read2.qname," with output ", output
                return -1
            else:
                tabWriter.writerow(output)
                return 1
        # CaseD: < 2 CpGs in read-pair
        elif (len(CpG_index1) + len(CpG_index2)) < 2:
            return 1
        else:
            print 'ERROR: Unexpected CpG count for read-pair ', read1.qname
            return -1

        # Case2: read-pair informative for OB strand, read-pair orientation is -/+ (read1/read2)
    elif read1.is_reverse and not read2.is_reverse:

        # If read2.tlen < 0, then the read-pair has a bad overlap
        if read2.tlen < 0:
            print 'WARNING: Ignoring read-pair ', read2.qname, ' due to bad overlap (=', read2.tlen, ')'
            return 0
        
        # start1 >= start2 by definition for a properly paired OB read
        start1 = read1.pos + 1 # read1.pos is 0-based in accordance with BAM file specificiations, regardless of whether the input file is SAM or BAM. I convert this to a 1-based position.
        start2 = read2.pos + 1 # read2.pos is 0-based in accordance with BAM file specificiations, regardless of whether the input file is SAM or BAM. I convert this to a 1-based position.
        strand = "-"
        
        # If the two reads overlap, i.e the rightmost position of read2 is greater than or equal to the leftmost position of read1, drop the offending positions from the rightmost (3') position of read2
        if read2.positions[-1] >= read1.positions[0]:
            ignore2 = len(read2XM) - int(read2.positions[-1] - read1.positions[0] + 1) # In this case, ignore2 is the number of bases to retain from the 5' end of read2.
            read2XM = read2XM[:ignore2]
        
        # Identify CpGs in read1 and read2
        CpG_index1 = [m.start() for m in re.finditer(CpG_pattern, read1XM)]
        CpG_index2 = [m.start() for m in re.finditer(CpG_pattern, read2XM)]
        
        # CaseA: A CpG in both reads of pair
        if (len(CpG_index1) > 0 and len(CpG_index2) > 0):
            CpGL = CpG_index2[0] # Leftmost CpG in read2
            CpGR = CpG_index1[-1] # Rightmost CpG in read1
            positionL = start2 + CpGL - 1 # No need to ignore2 since 3' bases of read2 are rightmost; -1 so as to point to C on the forward strand in the CpG 
            positionR = start1 + CpGR - 1 # -1 so as to point to C on the forward strand in the CpG
            output=[chrom, positionL, positionR, read2XM[CpGL], read1XM[CpGR], strand] # Output is left-to-right, thus read2XM appears before read1XM for OB read-pairs.
            if positionL > positionR:
                print "ERROR: Case2A posL > posR for read-pair ", read1.qname," with output ", output
                return -1
            else:
                tabWriter.writerow(output)
                return 1
        # CaseB: > 1 CpG in read1 but 0 CpGs in read2
        elif (len(CpG_index1) > 1 and len(CpG_index2) == 0):
            CpGL = CpG_index1[0] # Leftmost CpG in read1
            CpGR = CpG_index1[-1] # Rightmost CpG in read1
            positionL = start1 + CpGL - 1 # -1 so as to point to C on the forward strand in the CpG 
            positionR = start1 + CpGR - 1 # -1 so as to point to C on the forward strand in the CpG
            output=[chrom, positionL, positionR, read1XM[CpGL], read1XM[CpGR], strand]
            if positionL > positionR:
                print "ERROR: Case2B posL > posR for read ", read1.qname," with output ", output
                return -1
            else:
                tabWriter.writerow(output)
                return 1
        # CaseC: 0 CpGs in read1 but > 1 CpG in read2
        elif (len(CpG_index1) == 0 and len(CpG_index2) > 1):
            CpGL = CpG_index2[0] # Leftmost CpG in read2
            CpGR = CpG_index2[-1] # Rightmost CpG in read2
            positionL = start2 + CpGL - 1 # -1 so as to point to C on the forward strand in the CpG 
            positionR = start2 + CpGR - 1 # -1 so as to point to C on the forward strand in the CpG
            output=[chrom, positionL, positionR, read2XM[CpGL], read2XM[CpGR], strand]
            if positionL > positionR:
                print "ERROR: Case2C posL > posR for read ", read2.qname," with output ", output
                return -1
            else:
                tabWriter.writerow(output)
                return 1
        # CaseD: < 2 CpGs in read-pair
        elif (len(CpG_index1) + len(CpG_index2)) < 2:
            return 1
        else:
            print 'ERROR: Unexpected CpG count for read-pair', read1.qname
            return -1
    # Case3: read-pair in unexpected orientation
    else:
        print 'ERROR: Unexpected orientation of read-pair', read1.qname
        return -1


## tabWriter writes a tab-separated output file
tabWriter = csv.writer(OUT, delimiter='\t', quotechar=' ', quoting=csv.QUOTE_MINIMAL)

##############################################################################################################################################################################################

### The main program. Loops over the BAM file line-by-line (i.e. alignedRead-by-alignedRead) and extracts the XM information for each read or read-pair. ###
##############################################################################################################################################################################################
for read in BAM:
    if args.ignoreDuplicates and read.is_duplicate:
        continue
    else:
        if read.is_read1 and read.is_proper_pair:
            read1 = read
            continue
        elif read.is_read2 and read.is_proper_pair: 
            read2 = read
            # Check that read1 and read2 are aligned to the same chromosome and have identical read-names. If not, skip the read-pair.
            if read1.tid == read2.tid and read1.qname == read2.qname:
                PE_return = PE_SAM2MS(read1, read2)
                if PE_return == 0:
                    bad_overlaps_pairs += 1
            elif read1.tid != read2.tid:
                print "ERROR: Reads aligned to different chromosomes: ", read1.qname, BAM.getrname(read1.tid), read2.qname, BAM.getrname(read2.tid)
                continue
            elif read1.qname != read2.qname:
                print "ERROR: The name of read1 is not identical to to that of read2 for read-pair ", read1.qname, read2.qname
                break
        elif read.is_paired and not read.is_proper_pair:
            print "WARNING: Skipping read ", read.qname, " as it is part of an improper pair."
            continue
        elif not read.is_paired:  
            print "SE_SAM2MS(read)"
        #else:
         #   print "Read is neither a single-end read nor part of a paired-end read!", read.qname   
        
BAM.close()
OUT.close()
##############################################################################################################################################################################################
